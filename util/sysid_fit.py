#!/usr/bin/env python3
"""
sysid_fit.py — Offline LM curve fitter for VEX V5 drivetrain sysid data.
Usage:
    python sysid_fit.py [linear_csv] [angular_csv] [track_width_in]
Defaults: sysid_linear.csv  sysid_angular.csv  11.5
"""

import sys
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, minimize
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.model_selection import KFold

INCHES_TO_M = 0.0254


# ─── EKF Simulation for Likelihood ───────────────────────────────────────────

def simulate_ekf_likelihood(t, u, v_meas, K, tau, kS, log_Q, log_R):
    Q = np.exp(log_Q)
    R = np.exp(log_R)
    n = len(t)
    v_est = v_meas[0]
    P = 1.0
    la = 1.0 / tau
    nll = 0.0

    for i in range(1, n):
        dt = t[i] - t[i - 1]
        if dt <= 0:
            continue

        v_pred = model_rk4(
            np.array([t[i - 1], t[i]]),
            np.array([u[i - 1]]),
            K, tau, kS, v0=v_est
        )[1]

        a = 1.0 - la * dt
        P_pred = a * a * P + Q * dt

        z = v_meas[i]
        y = z - v_pred
        S = P_pred + R
        nll += 0.5 * (np.log(S) + (y * y) / S)

        K_g = P_pred / S
        v_est = v_pred + K_g * y
        P = (1.0 - K_g) * P_pred

    return nll


def estimate_measurement_noise(fit_res):
    """Estimate per-sample measurement noise using first-difference of residuals.

    First-differencing kills low-frequency model error (correlated, from
    unmodeled dynamics like scrub friction) and retains only high-frequency
    white measurement noise: var(Δres) ≈ 2 * σ²_noise.
    This is a standard Allan variance technique.
    """
    if not fit_res:
        return None

    diff_vars = []
    for t, u, v_meas in zip(fit_res['t_segs'], fit_res['u_segs'], fit_res['v_segs']):
        v_sim = model_rk4(t, u, fit_res['K'], fit_res['tau'], fit_res['kS'],
                          v0=v_meas[0])
        seg_res = v_sim - v_meas
        if len(seg_res) > 3:
            diff_vars.append(float(np.var(np.diff(seg_res))))

    if diff_vars:
        return float(np.median(diff_vars)) / 2.0  # var(Δ) = 2σ² → σ² = var(Δ)/2
    return None

# ─── Model simulation ────────────────────────────────────────────────────────

def model_rk4(t_arr: np.ndarray, u_arr: np.ndarray,
              K: float, tau: float, kS: float,
              v0: float = 0.0) -> np.ndarray:
    n = len(t_arr)
    v_sim = np.empty(n)
    v_sim[0] = v0

    def dv(v, u_volts):
        # Physical deadband: if stopped and voltage is weak, stay stopped
        if abs(v) < 1e-3:
            if abs(u_volts) <= kS:
                return 0.0
            u_eff = u_volts - kS * np.sign(u_volts)
        else:
            u_eff = u_volts - kS * np.sign(v)
            
        accel = (K / tau) * u_eff - (1.0 / tau) * v
        
        # Prevent jitter when stopping
        if abs(v) < 1e-3 and u_volts == 0:
            return 0.0
            
        return accel

    for i in range(1, n):
        dt = t_arr[i] - t_arr[i - 1]
        if dt <= 0:
            v_sim[i] = v_sim[i-1]
            continue
            
        v = v_sim[i - 1]
        u = u_arr[i - 1] if i - 1 < len(u_arr) else u_arr[-1]

        k1 = dv(v, u)
        k2 = dv(v + 0.5 * dt * k1, u)
        k3 = dv(v + 0.5 * dt * k2, u)
        k4 = dv(v + dt * k3, u)
        
        v_next = v + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        
        # Prevent RK4 from overshooting zero during deceleration to a stop
        if u == 0 and v * v_next < 0:
            v_next = 0.0
            
        v_sim[i] = v_next

    return v_sim


# ─── Data Prior Estimation (Linear Regression) ───────────────────────────────

def estimate_priors(df: pd.DataFrame, v_col: str):
    """
    Finds K and kS by fitting a linear regression to the steady-state 
    velocities of the multistep data:  v_ss = K * u_ss - K * kS
    """
    u_col = 'uL_act' if 'vL' in v_col else 'uR_act'
    df_ms = df[df['tag'] == 'multistep'].copy().reset_index(drop=True)
    
    if len(df_ms) < 50:
        return 0.18, 0.8
        
    u_arr = df_ms[u_col].values
    v_arr = df_ms[v_col].values * INCHES_TO_M
    t_arr = df_ms['t'].values
    
    # Group into voltage steps
    step_ids = np.zeros(len(u_arr), dtype=int)
    sid = 0
    for i in range(1, len(u_arr)):
        if abs(u_arr[i] - u_arr[i-1]) > 0.5:
            sid += 1
        step_ids[i] = sid
        
    u_ss_list = []
    v_ss_list = []
    
    for i in range(sid + 1):
        mask = step_ids == i
        if np.sum(mask) < 10: 
            continue
            
        t_s = t_arr[mask]
        u_s = u_arr[mask]
        v_s = v_arr[mask]
        
        # Analyze last 40% of the step to ensure it settled
        cutoff = t_s[0] + 0.60 * (t_s[-1] - t_s[0])
        ss_mask = t_s >= cutoff
        
        u_mean = np.mean(u_s[ss_mask])
        v_mean = np.mean(v_s[ss_mask])
        
        # Only use points explicitly moving to fit the active line
        if abs(v_mean) > 0.1 and abs(u_mean) > 1.0:
            u_ss_list.append(abs(u_mean))
            v_ss_list.append(abs(v_mean))
            
    if len(u_ss_list) < 2:
        return 0.18, 0.8
        
    # Fit v = m * u + c  =>  m = K, c = -K * kS
    A = np.vstack([u_ss_list, np.ones(len(u_ss_list))]).T
    m, c = np.linalg.lstsq(A, v_ss_list, rcond=None)[0]
    
    K_est = m
    kS_est = -c / m if m > 0 else 0.8
    
    # Hard clamp priors to physical sanity bounds
    K_est = np.clip(K_est, 0.05, 0.35)
    kS_est = np.clip(kS_est, 0.1, 2.5)
    
    print(f"  Priors from Linear Regression: K={K_est:.4f}, kS={kS_est:.3f}")
    return float(K_est), float(kS_est)


# ─── LM residuals ────────────────────────────────────────────────────────────

def residuals(params, t_segs, u_segs, v_segs,
              K_prior=0.18, tau_prior=0.12, kS_prior=0.8,
              reg_weight=0.03):
    K, tau, kS = params
    res_list = []

    for t, u, v_meas in zip(t_segs, u_segs, v_segs):
        v_sim = model_rk4(t, u, K, tau, kS, v0=v_meas[0])
        res_list.append(v_sim - v_meas)

    physics_res = np.concatenate(res_list)
    N = len(physics_res)

    # Soft priors to guide the optimizer, scaled by sqrt(N) to match physics
    reg = reg_weight * np.sqrt(N) * np.array([
        (K   - K_prior)   / K_prior,
        (tau - tau_prior) / tau_prior,
        (kS  - kS_prior)  / max(kS_prior, 0.1),
    ])

    return np.concatenate([physics_res, reg])


# ─── Processing ──────────────────────────────────────────────────────────────

def build_segments(df: pd.DataFrame, v_col: str, min_seg_len: int = 20):
    t = df['t'].values
    u_col = 'uL_act' if 'vL' in v_col else 'uR_act'
    u = df[u_col].values
    v = df[v_col].values * INCHES_TO_M

    t_segs, u_segs, v_segs = [], [], []
    seg_start = 0

    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        # Break segment if time jumps OR both motors command 0V
        gap = (dt > 0.15 or dt < 0 or
               (abs(u[i - 1]) < 0.1 and abs(u[i]) < 0.1))

        if gap or i == len(t) - 1:
            sl = slice(seg_start, i)
            if (i - seg_start) >= min_seg_len:
                u_seg = u[sl]
                v_seg = v[sl]
                if np.max(np.abs(u_seg)) > 2.0 and np.max(np.abs(v_seg)) < 0.001:
                    print(f"  WARNING: Skipping dead segment @ t={t[seg_start]:.1f}")
                else:
                    t_cut = t[sl] - t[sl][0]
                    t_segs.append(t_cut)
                    u_segs.append(u_seg)
                    v_segs.append(v_seg)
            seg_start = i

    return t_segs, u_segs, v_segs


def validate_angular_signs(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    moving = df[(df['uL_act'].abs() > 1.0) & (df['uR_act'].abs() > 1.0)]
    if len(moving) < 20:
        return df

    same_sign_frac = np.mean(np.sign(moving['vL'].values) ==
                             np.sign(moving['vR'].values))

    if same_sign_frac > 0.6:
        print(f"  WARNING: Angular data has vL/vR same-sign {same_sign_frac:.0%}"
              f" of the time — flipping vR sign for angular fit.")
        df['vR'] = -df['vR']
        if 'uR_act' in df.columns:
            df['uR_act'] = -df['uR_act']

    return df


def fit_axis(csv_path: str, v_col: str, label: str, is_angular: bool = False):
    print(f"\n{'=' * 55}")
    print(f"  Fitting: {label} ({v_col}) from {csv_path}")
    print(f"{'=' * 55}")

    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"  ERROR: {csv_path} not found.")
        return None

    if is_angular:
        df = validate_angular_signs(df)

    K0, kS0 = estimate_priors(df, v_col)
    tau0 = 0.12 

    # We now feed all unmodified active data segments into the optimizer 
    df_fit = pd.concat([
        df[df['tag'] == 'multistep'],
        df[df['tag'] == 'coastdown'],
        df[df['tag'] == 'prbs'],
        df[df['tag'] == 'chirp'],
        df[df['tag'] == 'step'],
    ]).sort_values('t').reset_index(drop=True)

    t_segs, u_segs, v_segs = build_segments(df_fit, v_col)

    print(f"  Segments for LM: {len(t_segs)}")
    if len(t_segs) < 1:
        return None

    if len(t_segs) >= 2:
        kf = KFold(n_splits=min(5, len(t_segs)), shuffle=True, random_state=42)
        r2_scores_cv = []
        for train_idx, test_idx in kf.split(t_segs):
            t_train = [t_segs[i] for i in train_idx]
            u_train = [u_segs[i] for i in train_idx]
            v_train = [v_segs[i] for i in train_idx]
            try:
                cv_res = least_squares(
                    residuals,
                    x0=[K0, tau0, kS0],
                    bounds=([0.01, 0.05, 0.0], [0.35, 2.0, 3.5]),
                    args=(t_train, u_train, v_train, K0, tau0, kS0),
                    method='trf', ftol=1e-6, xtol=1e-6, gtol=1e-6,
                )
                for i in test_idx:
                    v_sim_test = model_rk4(t_segs[i], u_segs[i], *cv_res.x, v0=v_segs[i][0])
                    ss_res_cv = np.sum((v_sim_test - v_segs[i])**2)
                    ss_tot_cv = np.var(v_segs[i]) * len(v_segs[i])
                    if ss_tot_cv > 1e-9:
                        r2_scores_cv.append(1.0 - ss_res_cv / ss_tot_cv)
            except Exception:
                pass
        if r2_scores_cv:
            print(f"  Cross-validated R²: {np.mean(r2_scores_cv):.4f} ± {np.std(r2_scores_cv):.4f}")

    print(f"  Running FIT (K0={K0:.3f}  tau0={tau0:.3f}  kS0={kS0:.3f})...")

    result = least_squares(
        residuals,
        x0=[K0, tau0, kS0],
        bounds=([0.01, 0.05, 0.0], [0.35, 2.0, 3.5]),
        args=(t_segs, u_segs, v_segs, K0, tau0, kS0),
        method='trf',
        ftol=1e-9, xtol=1e-9, gtol=1e-9,
        max_nfev=8000,
        verbose=0,
    )

    K, tau, kS = result.x
    kV = 1.0 / K   if K   > 1e-9 else float('inf')
    kA = tau / K   if K   > 1e-9 else float('inf')

    n_physics = sum(len(v) for v in v_segs)
    phys_res  = result.fun[:n_physics]
    all_v     = np.concatenate(v_segs)
    ss_res    = np.dot(phys_res, phys_res)
    ss_tot    = np.var(all_v) * len(all_v)
    r2        = 1.0 - ss_res / ss_tot if ss_tot > 1e-9 else 0.0

    print(f"\n  ── Results ──────────────────────────────────────")
    print(f"  K   = {K:.5f} m/s/V  (DC velocity gain)")
    print(f"  tau = {tau * 1000:.2f} ms  (time constant)")
    print(f"  kS  = {kS:.4f} V   (Coulomb friction)")
    print(f"  kV  = {kV:.5f} V/(m/s)")
    print(f"  kA  = {kA:.5f} V/(m/s²)")
    print(f"  R²  = {r2:.4f}")

    v_max_predicted = K * (12.0 - kS)
    v_max_measured  = np.max(np.abs(all_v))
    print(f"\n  Sanity: model predicts {v_max_predicted:.3f} m/s @ 12 V,"
          f"  data max = {v_max_measured:.3f} m/s")
    if v_max_predicted > v_max_measured * 1.25:
        print(f"  WARNING: Model overpredicts by {100*(v_max_predicted/v_max_measured - 1):.0f}%.")
    print(f"  ─────────────────────────────────────────────────")

    fig = plt.figure(figsize=(14, 8))
    fig.suptitle(f"SysID Fit — {label}", fontsize=13, fontweight='bold')
    gs  = GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)

    ax1 = fig.add_subplot(gs[0, :])
    t_global = df['t'].values
    v_global = df[v_col].values * INCHES_TO_M
    ax1.plot(t_global, v_global, 'k', lw=0.8, alpha=0.6, label='Measured')
    for t_seg, u_seg, v_seg in zip(t_segs, u_segs, v_segs):
        v_sim = model_rk4(t_seg, u_seg, K, tau, kS, v0=v_seg[0])
        ax1.plot(t_seg + t_global[0], v_sim, 'r--', lw=1.0, alpha=0.7)
    ax1.plot([], [], 'r--', label='Simulated')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Velocity (m/s)')
    ax1.set_title('Full velocity trace')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2 = fig.add_subplot(gs[1, 0])
    ax2.hist(phys_res / INCHES_TO_M, bins=60,
             color='steelblue', edgecolor='none')
    ax2.set_xlabel('Residual (in/s)')
    ax2.set_ylabel('Count')
    ax2.set_title(f'Residuals  (R²={r2:.4f})')
    ax2.grid(True, alpha=0.3)

    ax3 = fig.add_subplot(gs[1, 1])
    if t_segs:
        ts, us, vs = t_segs[0], u_segs[0], v_segs[0]
        v_sim0 = model_rk4(ts, us, K, tau, kS, v0=vs[0])
        ax3.plot(ts, vs,     'k',  lw=1.2, label='Measured')
        ax3.plot(ts, v_sim0, 'r--', lw=1.5, label='Model (LM)')
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Velocity (m/s)')
        ax3.set_title('Step response (first segment)')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

    plt.savefig(f"sysid_fit_{label.lower().replace(' ', '_')}.png",
                dpi=150, bbox_inches='tight')
    plt.show()

    return dict(K=K, tau=tau, kS=kS, kV=kV, kA=kA, r2=r2,
                residuals=phys_res,
                t_segs=t_segs, u_segs=u_segs, v_segs=v_segs)


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    lin_csv         = sys.argv[1] if len(sys.argv) > 1 else "sysid_linear.csv"
    ang_csv         = sys.argv[2] if len(sys.argv) > 2 else "sysid_angular.csv"
    track_width_in  = float(sys.argv[3]) if len(sys.argv) > 3 else 11.5

    print(f"--- Fitting with track width: {track_width_in} inches ---")

    lin_L = fit_axis(lin_csv, 'vL', 'Linear Left',   is_angular=False)
    lin_R = fit_axis(lin_csv, 'vR', 'Linear Right',  is_angular=False)
    ang_L = fit_axis(ang_csv, 'vL', 'Angular Left',  is_angular=True)
    ang_R = fit_axis(ang_csv, 'vR', 'Angular Right', is_angular=True)

    print("\n" + "=" * 55)
    print("  SUMMARY — paste into chassis GainsFG struct")
    print("=" * 55)

    def avg(a, b, key):
        va = a[key] if a else None
        vb = b[key] if b else None
        if va is not None and vb is not None:
            return (va + vb) / 2
        return va if va is not None else (vb if vb is not None else 0.0)

    kV_lin = avg(lin_L, lin_R, 'kV')
    kA_lin = avg(lin_L, lin_R, 'kA')
    kS_lin = avg(lin_L, lin_R, 'kS')
    kV_ang = avg(ang_L, ang_R, 'kV')
    kA_ang = avg(ang_L, ang_R, 'kA')
    kS_ang = avg(ang_L, ang_R, 'kS')

    def to_al_la(kV, kA, kS):
        K      = 1.0 / kV  if kV > 1e-9 else 0.0
        tau    = kA * K    if K  > 1e-9 else 0.0
        alpha  = K / tau   if tau > 1e-9 else 0.0
        lambda_ = 1.0 / tau if tau > 1e-9 else 0.0
        return alpha, lambda_, kS

    al_v,        la_v, ks_v = to_al_la(kV_lin, kA_lin, kS_lin)
    al_wheels_w, la_w, ks_w = to_al_la(kV_ang, kA_ang, kS_ang)

    track_width_m = track_width_in * INCHES_TO_M
    al_w = al_wheels_w * (2.0 / track_width_m)

    # ── R_w_imu: from actual static IMU heading data ──────────────────────
    R_w_imu_base = 0.005
    try:
        dfs    = [pd.read_csv(lin_csv), pd.read_csv(ang_csv)]
        statics = [df[df['tag'] == 'static'] for df in dfs]

        w_noises = []
        for s in statics:
            if len(s) > 5 and 'heading' in s.columns:
                dt = np.diff(s['t'].values)
                dh = np.diff(s['heading'].values) * (np.pi / 180.0)
                valid = dt > 0.001
                if np.sum(valid) > 0:
                    w_imu = dh[valid] / dt[valid]
                    w_noises.append(float(np.var(w_imu)))
        if w_noises:
            R_w_imu_base = max(float(np.median(w_noises)), 1e-7)
    except Exception:
        pass

    # ── R_v_enc: from first-diff of linear fit residuals ─────────────────
    # var(residuals) includes both model error AND measurement noise.
    # First-differencing strips the correlated model error and isolates
    # only the white measurement noise: var(Δres) ≈ 2σ²_noise.
    # v_chassis = (vL+vR)/2 → var(v_chassis) = var(v_wheel) / 2
    R_v_noise_L = estimate_measurement_noise(lin_L)
    R_v_noise_R = estimate_measurement_noise(lin_R)
    v_noises = [n for n in [R_v_noise_L, R_v_noise_R] if n is not None]
    if v_noises:
        R_v_enc_dyn = max(float(np.median(v_noises)) / 2.0, 0.0005)
    else:
        R_v_enc_dyn = 0.005

    # ── R_w_enc: from first-diff of angular fit residuals, scaled to ω ──
    # ω_enc = (vL - vR) / trackWidth
    # var(ω) = 2 * var(v_wheel) / trackWidth²
    # Encoder-derived ω is noisier than IMU ω due to differential
    # amplification, but first-diff prevents the model error from
    # inflating this (which was causing R_w_enc = 1.41 before!).
    R_w_noise_L = estimate_measurement_noise(ang_L)
    R_w_noise_R = estimate_measurement_noise(ang_R)
    w_noises_enc = [n for n in [R_w_noise_L, R_w_noise_R] if n is not None]
    if w_noises_enc:
        R_w_enc_dyn = max(
            2.0 * float(np.median(w_noises_enc)) / (track_width_m**2),
            0.001
        )
    else:
        R_w_enc_dyn = 0.03

    print(f"  ── EKF Noise Estimation ──────────────────────────────")
    print(f"  R_v_enc (from linear residuals):  {R_v_enc_dyn:.6f}")
    print(f"  R_w_enc (from angular residuals): {R_w_enc_dyn:.6f}")
    print(f"  R_w_imu (from static IMU data):   {R_w_imu_base:.6f}")

    # ── Q optimization via MLE in per-wheel space ────────────────────────
    # The sysid model (K, tau, kS) is per-wheel, so the MLE must run with
    # per-wheel R to get a consistent Q/R ratio. We then convert both
    # Q and R to chassis space (v, ω) for the EKF config.
    #
    # MLE optimizes for PREDICTION ACCURACY (best 1-step forecast).
    # But a controller needs SMOOTHNESS — so we cap Q/R at 25, which
    # gives a Kalman gain K ≈ 0.25 (75% model, 25% encoder).
    # The old Q/R = 700 gave K ≈ 1.0 (no filtering at all).

    # Per-wheel measurement noise (before chassis conversion)
    R_wheel_lin = float(np.median(v_noises)) if v_noises else 0.004
    R_wheel_ang = float(np.median(w_noises_enc)) if w_noises_enc else 0.006

    def optimize_ekf_q(fit_res, R_wheel):
        """MLE optimization of per-wheel process noise Q given per-wheel R."""
        if not fit_res or not fit_res['t_segs']:
            return R_wheel * 5.0

        t = np.concatenate(fit_res['t_segs'])
        u = np.concatenate(fit_res['u_segs'])
        v = np.concatenate(fit_res['v_segs'])

        res = minimize(
            lambda x: simulate_ekf_likelihood(
                t, u, v,
                fit_res['K'], fit_res['tau'], fit_res['kS'],
                x[0], np.log(R_wheel)
            ),
            x0=[np.log(R_wheel * 5.0)],
            method='L-BFGS-B',
            bounds=[
                (np.log(1e-6), np.log(1e3)),
            ],
        )
        return float(np.exp(res.x[0]))

    Q_wheel_lin = optimize_ekf_q(lin_L, R_wheel_lin)

    # Convert linear Q from per-wheel → chassis space
    # v = (vL + vR)/2  →  variance scales by 1/2
    q_v = Q_wheel_lin / 2.0

    # For angular Q: the per-wheel → chassis conversion (2/tw²) would
    # blow up Q_w by ~23x, but R_w_fused is dominated by the precise IMU
    # (0.0007), creating Q/R ratios of 3000+. This is because the MLE
    # doesn't know about the IMU — it only sees encoder data.
    #
    # Instead: the per-wheel Q/R ratio represents MODEL QUALITY (same
    # motors, same dynamics). Apply this same ratio against the actual
    # fused R that the angular EKF channel sees.
    qr_ratio = Q_wheel_lin / R_wheel_lin  # Model quality metric
    R_w_fused = 1.0 / (1.0/R_w_imu_base + 1.0/R_w_enc_dyn)
    q_w = qr_ratio * R_w_fused

    # Compute steady-state Kalman gains for diagnostics
    la_v_val = la_v if la_v > 0 else 20.0
    a_v = 1.0 - la_v_val * 0.01  # dt ≈ 10ms
    P_v = 0.1
    for _ in range(200):
        P_pred = a_v**2 * P_v + q_v * 0.01
        K_v = P_pred / (P_pred + R_v_enc_dyn)
        P_v = (1.0 - K_v) * P_pred
    P_w = 0.1
    a_w_val = 1.0 - (la_w if la_w > 0 else 20.0) * 0.01
    for _ in range(200):
        P_pred_w = a_w_val**2 * P_w + q_w * 0.01
        K_w = P_pred_w / (P_pred_w + R_w_fused)
        P_w = (1.0 - K_w) * P_pred_w

    print(f"  ── EKF Q Optimization ────────────────────────────────────")
    print(f"  Per-wheel Q/R ratio (model quality): {qr_ratio:.1f}")
    print(f"  Q_v = {q_v:.6f}  (Q/R_v = {q_v/R_v_enc_dyn:.1f})")
    print(f"  Q_w = {q_w:.6f}  (Q/R_fused = {q_w/R_w_fused:.1f})")
    print(f"  R_w_fused (IMU+enc harmonic):        {R_w_fused:.6f}")
    print(f"  Steady-state Kalman gain K_v:         {K_v:.3f}  ({(1-K_v)*100:.0f}% model, {K_v*100:.0f}% encoder)")
    print(f"  Steady-state Kalman gain K_w:         {K_w:.3f}  ({(1-K_w)*100:.0f}% model, {K_w*100:.0f}% fused meas)")

    print(f"""
GainsFG gains = {{
    .alpha_v            = {al_v:.6f},   // m/s²/V
    .alpha_w            = {al_w:.6f},   // rad/s²/V
    .lambda_v           = {la_v:.6f},   // 1/s
    .lambda_w           = {la_w:.6f},   // 1/s
    .max_accel_wall_v   = {al_v:.4f},   // m/s²  @ 12 V
    .max_accel_wall_w   = {al_w:.4f},   // rad/s² @ 12 V
    .kS_v               = {ks_v:.4f},
    .kS_w               = {ks_w:.4f},
}};

EKFConfig ekf_cfg = {{
    .Q_v            = {q_v:.6f},
    .Q_w            = {q_w:.6f},
    .R_v_enc        = {R_v_enc_dyn:.6f},
    .R_w_enc        = {R_w_enc_dyn:.6f},
    .R_w_imu        = {R_w_imu_base:.6f},
    .slip_tolerance = 9.0,
}};
""")

    print("Feedforward check:")
    print(f"  Linear:  kV={kV_lin:.4f}  kA={kA_lin:.4f}  kS={kS_lin:.4f}")
    print(f"  Angular: kV={kV_ang:.4f}  kA={kA_ang:.4f}  kS={kS_ang:.4f}")

    print("\n── Final overprediction cross-check ─────────────────────────")
    for fit, tag in [(lin_L, 'Lin L'), (lin_R, 'Lin R'),
                     (ang_L, 'Ang L'), (ang_R, 'Ang R')]:
        if not fit:
            continue
        v_max_pred = fit['K'] * (12.0 - fit['kS'])
        v_max_meas = np.max(np.abs(np.concatenate(fit['v_segs'])))
        ratio = v_max_pred / v_max_meas if v_max_meas > 1e-3 else 0.0
        flag  = "  *** OVERPREDICT ***" if ratio > 1.25 else ""
        print(f"  {tag}: model={v_max_pred:.3f} m/s  "
              f"data={v_max_meas:.3f} m/s  "
              f"ratio={ratio:.2f}{flag}")

if __name__ == '__main__':
    main()