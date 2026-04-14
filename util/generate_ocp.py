from acados_template import AcadosOcp, AcadosModel, AcadosOcpSolver
from casadi import SX, vertcat, cos, sin, tanh, fmax, if_else, logic_and, sqrt, atan2
import numpy as np

"""
Simplified Path-Following MPC for Differential Drive

States: [x, y, theta, vy, s] - position, heading, lateral velocity, path progress
Controls: [v, w, vs] - forward velocity, angular velocity, path progress rate

Key features:
1. Continuous heading cost using (cos(th)-cos(ref))² + (sin(th)-sin(ref))²
   - No wrap discontinuity at ±180°!
   - Equivalent to 2(1 - cos(error)) ≈ error² for small errors
2. Cross-track error for speed vs accuracy tuning
3. Keeps s as state for natural speed adaptation in curves

The path is provided as a waypoint buffer with Hermite spline interpolation.
"""

def build_solver_local(
    solver_name="mpcc",
    dt_list=None,
    # Tracking weights
    q_pos=10.0,      # Position error weight (Euclidean distance)
    q_heading=5.0,   # Heading error weight (continuous formula)
    q_cross=8.0,     # Cross-track error weight (sideways from path)
    q_vy=4.0,        # Lateral velocity penalty
    # Progress reward (encourages forward motion along path)
    q_progress=2.0,
    # Control rate penalties (smooth control changes)
    w_dv=5.0,
    w_dw=5.0,
    w_dvs=3.0,
    # Velocity tracking
    w_v=6.0,
    v_ref=0.8,
    # Terminal weights
    q_pos_N=20.0,
    q_heading_N=10.0,
    q_cross_N=16.0,
    q_vy_N=8.0,
    # Control bounds
    v_min=0.0,
    v_max=1.5,
    w_min=-2.5,
    w_max=2.5,
    vs_max=2.0,
    # Physical parameters
    fk=0.3,          # Friction coefficient for drift model
    d_intake=0.25,   # Distance from center to intake point
    s_trust=0.3,     # Trust region half-width for progress
    # Curvature-based adaptation
    kappa_threshold=4.0,
):
    if dt_list is None:
        dt_list = [0.01]*10 + [0.02]*15

    N = len(dt_list)
    T = float(np.sum(dt_list))

    # ============================================================
    # STATE DEFINITION: [x, y, theta, vy, s]
    # ============================================================
    x_ = SX.sym("x")
    y_ = SX.sym("y")
    th = SX.sym("theta")  # Compass heading: 0=North, CW positive
    vy = SX.sym("vy")     # Lateral velocity (drift)
    s = SX.sym("s")       # Path progress (arc length)
    x_state = vertcat(x_, y_, th, vy, s)

    # ============================================================
    # CONTROL DEFINITION: [v, w, vs]
    # ============================================================
    v = SX.sym("v")   # Forward velocity (m/s)
    w = SX.sym("w")   # Angular velocity (rad/s), CW positive
    vs = SX.sym("vs") # Path progress rate (m/s)
    u_control = vertcat(v, w, vs)

    # ============================================================
    # PARAMETERS (78 total - added q_cross)
    # ============================================================
    # p[0]: dt
    # p[1]: fk (friction)
    # p[2]: d_intake
    # p[3-5]: prev controls [v, w, vs]
    # p[6]: s_center (for trust region)
    # p[7-13]: weights [q_pos, q_heading, q_cross, q_vy, q_progress, w_dv, w_dw]
    # p[14-19]: more weights [w_dvs, w_v, v_ref, q_pos_N, q_heading_N, q_cross_N]
    # p[20]: q_vy_N
    # p[21]: kappa (curvature at this node)
    # p[22]: kappa_threshold
    # p[23]: s_trust
    # p[24-31]: UNUSED (reserved)
    # p[32-76]: Waypoint buffer (9 waypoints * 5 components)

    p = SX.sym("p", 77)

    # Extract parameters
    dt_p = p[0]
    fk_p = p[1]
    d_intake_p = p[2]
    v_prev = p[3]
    w_prev = p[4]
    vs_prev = p[5]
    s_center = p[6]

    q_pos_p = p[7]
    q_heading_p = p[8]
    q_cross_p = p[9]      # NEW: cross-track weight
    q_vy_p = p[10]
    q_progress_p = p[11]
    w_dv_p = p[12]
    w_dw_p = p[13]
    w_dvs_p = p[14]
    w_v_p = p[15]
    v_ref_p = p[16]
    q_pos_N_p = p[17]
    q_heading_N_p = p[18]
    q_cross_N_p = p[19]   # NEW: terminal cross-track weight
    q_vy_N_p = p[20]
    kappa = p[21]
    kappa_threshold_p = p[22]
    s_trust_p = p[23]

    # ============================================================
    # DYNAMICS (Compass heading: theta=0 North, CW positive)
    # ============================================================
    # Midpoint integration for theta (more accurate for rotations)
    th_mid = th + 0.5 * w * dt_p

    # Position update: robot moves in direction (theta)
    # In compass: sin(th) = dx/ds, cos(th) = dy/ds
    dx = (v * sin(th_mid) + vy * cos(th_mid)) * dt_p
    dy = (v * cos(th_mid) - vy * sin(th_mid)) * dt_p
    dth = w * dt_p

    # Lateral velocity (drift) model - decays with friction
    # Centripetal force creates drift during turns
    f_k_dt = fk_p * dt_p
    centripetal = v * w * dt_p
    cent_clamped = f_k_dt * tanh(centripetal / (f_k_dt + 1e-6))
    vy_sign = tanh(100.0 * vy)  # Smooth sign function
    vy_next = vy - cent_clamped - vy_sign * f_k_dt

    # Path progress update
    s_next = s + vs * dt_p

    # Assemble state transition
    x_next = vertcat(x_ + dx, y_ + dy, th + dth, vy_next, s_next)

    # ============================================================
    # PATH REFERENCE INTERPOLATION (Hermite spline)
    # ============================================================
    # Waypoint buffer starts at p[32], each waypoint has 5 components
    W_OFFSET = 32
    W_SIZE = 5
    NUM_WAYPOINTS = 9

    # Default values (will be overwritten by segment search)
    xr = 0.0
    yr = 0.0
    ref_th = 0.0

    # Search through 8 segments to find which one contains current s
    for i in range(NUM_WAYPOINTS - 1):
        s_i = p[W_OFFSET + i * W_SIZE + 4]
        s_i1 = p[W_OFFSET + (i + 1) * W_SIZE + 4]

        # Cubic Hermite interpolation
        L = s_i1 - s_i
        t = (s - s_i) / fmax(L, 0.001)

        # Position interpolation
        p0_x = p[W_OFFSET + i * W_SIZE + 0]
        p0_y = p[W_OFFSET + i * W_SIZE + 1]
        v0_x = p[W_OFFSET + i * W_SIZE + 2]
        v0_y = p[W_OFFSET + i * W_SIZE + 3]
        p1_x = p[W_OFFSET + (i + 1) * W_SIZE + 0]
        p1_y = p[W_OFFSET + (i + 1) * W_SIZE + 1]
        v1_x = p[W_OFFSET + (i + 1) * W_SIZE + 2]
        v1_y = p[W_OFFSET + (i + 1) * W_SIZE + 3]

        # Hermite basis functions
        h00 = 2*t**3 - 3*t**2 + 1
        h10 = t**3 - 2*t**2 + t
        h01 = -2*t**3 + 3*t**2
        h11 = t**3 - t**2

        # Interpolated position
        x_interp = h00*p0_x + h10*L*v0_x + h01*p1_x + h11*L*v1_x
        y_interp = h00*p0_y + h10*L*v0_y + h01*p1_y + h11*L*v1_y

        # Derivative of Hermite basis (for tangent/heading)
        hp00 = 6*t**2 - 6*t
        hp10 = 3*t**2 - 4*t + 1
        hp01 = -6*t**2 + 6*t
        hp11 = 3*t**2 - 2*t

        # Tangent direction
        tx = hp00*p0_x + hp10*L*v0_x + hp01*p1_x + hp11*L*v1_x
        ty = hp00*p0_y + hp10*L*v0_y + hp01*p1_y + hp11*L*v1_y

        # Heading from tangent (compass convention: atan2(dx, dy))
        th_interp = atan2(tx / fmax(L, 0.001), ty / fmax(L, 0.001))

        # Use this segment if s is within it
        in_segment = logic_and(s >= s_i, s < s_i1)
        xr = if_else(in_segment, x_interp, xr)
        yr = if_else(in_segment, y_interp, yr)
        ref_th = if_else(in_segment, th_interp, ref_th)

    # Clamp to last waypoint if past the end
    final_s = p[W_OFFSET + (NUM_WAYPOINTS - 1) * W_SIZE + 4]
    xr = if_else(s >= final_s, p[W_OFFSET + (NUM_WAYPOINTS - 1) * W_SIZE + 0], xr)
    yr = if_else(s >= final_s, p[W_OFFSET + (NUM_WAYPOINTS - 1) * W_SIZE + 1], yr)
    ref_th = if_else(s >= final_s, p[W_OFFSET + (NUM_WAYPOINTS - 1) * W_SIZE + 2], ref_th)

    # ============================================================
    # ERROR COMPUTATION - CONTINUOUS ACROSS ±180°
    # ============================================================
    # Use intake point for error calculation (where the "action" happens)
    xi = x_ + d_intake_p * sin(th)
    yi = y_ + d_intake_p * cos(th)

    # 1. Euclidean distance to reference (linear cost for smooth degradation)
    pos_error = sqrt((xi - xr)**2 + (yi - yr)**2 + 1e-6)  # Added eps to avoid gradient singularity

    # 2. CONTINUOUS HEADING COST
    # (cos(th) - cos(ref))² + (sin(th) - sin(ref))² = 2(1 - cos(th - ref))
    # This is continuous everywhere (no ±180° discontinuity!)
    # For small errors: ≈ (th - ref)² since 1 - cos(x) ≈ x²/2
    heading_cost = (cos(th) - cos(ref_th))**2 + (sin(th) - sin(ref_th))**2

    # 3. CROSS-TRACK ERROR (signed distance perpendicular to path)
    # Path tangent direction (compass: sin=dx, cos=dy)
    path_tangent_x = sin(ref_th)
    path_tangent_y = cos(ref_th)
    # Path normal (perpendicular, pointing LEFT of forward direction)
    path_normal_x = cos(ref_th)
    path_normal_y = -sin(ref_th)
    # Signed cross-track error: positive = left of path, negative = right
    cross_track = (xi - xr) * path_normal_x + (yi - yr) * path_normal_y

    # 4. Along-track error (how far ahead/behind on path)
    along_track = (xi - xr) * path_tangent_x + (yi - yr) * path_tangent_y

    # ============================================================
    # CURVATURE-BASED ADAPTATION
    # ============================================================
    # Reduce target speed in high-curvature sections
    v_safe = v_ref_p / (1.0 + 0.5 * kappa / fmax(kappa_threshold_p, 0.5))

    # Increase cross-track weight in curves (need tighter path following)
    q_cross_adaptive = q_cross_p * (1.0 + kappa / fmax(kappa_threshold_p, 0.5))

    # ============================================================
    # STAGE COST
    # ============================================================
    cost_stage = (
        # Position tracking - Euclidean distance (linear, less aggressive)
        q_pos_p * pos_error

        # Heading tracking - CONTINUOUS cost (no wrap issues!)
        + 0.5 * q_heading_p * heading_cost

        # Cross-track error - penalize being off the path sideways
        # Higher q_cross = more accurate path following
        # Lower q_cross = faster but cuts corners
        + 0.5 * q_cross_adaptive * cross_track**2

        # Lateral velocity penalty (discourage drift)
        + 0.5 * q_vy_p * vy**2

        # Progress reward (linear, encourages forward motion)
        - q_progress_p * vs

        # Control rate penalties (smooth transitions)
        + 0.5 * w_dv_p * (v - v_prev)**2
        + 0.5 * w_dw_p * (w - w_prev)**2
        + 0.5 * w_dvs_p * (vs - vs_prev)**2

        # Velocity tracking (soft constraint on speed)
        + 0.5 * w_v_p * (v - v_safe)**2
    )

    # ============================================================
    # TRUST REGION (Soft constraint on progress deviation)
    # ============================================================
    s_upper = s_center + s_trust_p
    s_lower = s_center - s_trust_p
    s_viol = if_else(s > s_upper, s - s_upper,
                     if_else(s < s_lower, s_lower - s, 0.0))
    cost_stage += 50.0 * s_viol**2

    # ============================================================
    # TERMINAL COST
    # ============================================================
    cost_e = (
        q_pos_N_p * pos_error
        + 0.5 * q_heading_N_p * heading_cost
        + 0.5 * q_cross_N_p * cross_track**2
        + 0.5 * q_vy_N_p * vy**2
    )

    # ============================================================
    # BUILD MODEL
    # ============================================================
    model = AcadosModel()
    model.name = solver_name
    model.x = x_state
    model.u = u_control
    model.p = p
    model.disc_dyn_expr = x_next
    model.cost_expr_ext_cost = cost_stage
    model.cost_expr_ext_cost_e = cost_e

    # ============================================================
    # BUILD OCP
    # ============================================================
    ocp = AcadosOcp()
    ocp.model = model
    ocp.solver_options.N_horizon = N
    ocp.cost.cost_type = "EXTERNAL"
    ocp.cost.cost_type_e = "EXTERNAL"

    # Control bounds
    ocp.constraints.idxbu = np.array([0, 1, 2], dtype=int)
    ocp.constraints.lbu = np.array([v_min, w_min, 0.0])
    ocp.constraints.ubu = np.array([v_max, w_max, vs_max])
    ocp.constraints.x0 = np.zeros(5)

    # Default parameters
    p_def = np.zeros(77)
    p_def[0] = dt_list[0]
    p_def[1] = fk
    p_def[2] = d_intake
    p_def[7] = q_pos
    p_def[8] = q_heading
    p_def[9] = q_cross
    p_def[10] = q_vy
    p_def[11] = q_progress
    p_def[12] = w_dv
    p_def[13] = w_dw
    p_def[14] = w_dvs
    p_def[15] = w_v
    p_def[16] = v_ref
    p_def[17] = q_pos_N
    p_def[18] = q_heading_N
    p_def[19] = q_cross_N
    p_def[20] = q_vy_N
    p_def[22] = kappa_threshold
    p_def[23] = s_trust
    ocp.parameter_values = p_def

    # Solver options
    ocp.solver_options.time_steps = np.array(dt_list)
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    ocp.solver_options.integrator_type = "DISCRETE"
    ocp.solver_options.nlp_solver_type = "SQP_RTI"
    ocp.solver_options.regularize_method = "PROJECT"
    ocp.solver_options.qp_solver_cond_N = 5
    ocp.solver_options.qp_solver_iter_max = 50
    ocp.solver_options.qp_solver_warm_start = 2
    ocp.solver_options.ext_fun_compile_flags = "-O3"
    ocp.solver_options.tf = T
    ocp.solver_options.print_level = 0

    solver = AcadosOcpSolver(ocp, json_file=f"{solver_name}.json")

    # Initialize with reasonable guesses
    x_init = np.zeros(5)
    u_init = np.array([0.5, 0.0, 0.5])
    p_init = p_def.copy()

    for i in range(N + 1):
        solver.set(i, "x", x_init)
    for i in range(N):
        solver.set(i, "u", u_init)
        solver.set(i, "p", p_init)

    return solver, {
        "N": N,
        "dyn_params": {"fk": fk, "d_intake": d_intake},
        "solver_name": solver_name
    }


if __name__ == "__main__":
    solver, defaults = build_solver_local()
    print(f"Path MPC Solver Built: {defaults['solver_name']}")
    print(f"Horizon: {defaults['N']} steps")
