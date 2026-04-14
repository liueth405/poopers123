import argparse
import pandas as pd
import numpy as np
from scipy.optimize import minimize_scalar

def anglewrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi

def simulate(fk, data):
    n = len(data)
    time = data['time'].values / 1000.0 # seconds
    vl = data['vl'].values
    vr = data['vr'].values
    heading = data['curr_heading'].values
    
    dt_arr = np.diff(time)
    dt_arr = np.insert(dt_arr, 0, 0.02)
    
    x_sim = np.zeros(n)
    y_sim = np.zeros(n)
    theta_sim = np.zeros(n)
    
    x_sim[0] = data['pf_x'].iloc[0]
    y_sim[0] = data['pf_y'].iloc[0]
    theta_sim[0] = data['pf_theta'].iloc[0]
    
    lv = 0.0
    max_drift_speed = 41.0 # Derived from Particle Filter configs
    
    for i in range(1, n):
        dt = dt_arr[i]
        if dt <= 0:
            dt = 0.02
            
        dist = (vl[i] + vr[i]) * 0.5 * dt
        v_fwd = dist / dt if dt > 0 else 0.0
        
        dh = heading[i] - heading[i-1]
        dh = anglewrap(dh)
        
        f_k_dt = fk * dt
        
        # Calculate friction application given odom.hpp model
        raw_kick = v_fwd * dh
        kick_clamped = np.clip(raw_kick, -f_k_dt, f_k_dt)
        centrip_kick = kick_clamped
        
        v_pre_fric = lv - centrip_kick
        
        # Apply kinetic friction
        v_abs_pre = abs(v_pre_fric)
        v_reduced = v_abs_pre - f_k_dt
        v_clipped = max(0.0, v_reduced)
        
        lv_raw = v_clipped if v_pre_fric > 0 else -v_clipped
        
        max_vy_mag = abs(v_fwd)
        
        # Clamp bounds
        if abs(lv_raw) > max_vy_mag:
            lv_final = max_vy_mag if lv_raw > 0 else -max_vy_mag
        else:
            lv_final = lv_raw
            
        # Deadband
        if abs(lv_final) < 0.2:
            lv_new = 0.0
        else:
            lv_new = lv_final
            
        # Global clamp
        lv_new = np.clip(lv_new, -max_drift_speed, max_drift_speed)
        
        # Accumulate slip vector laterally
        deltaX = ((lv + lv_new) * 0.5) * dt
        
        # Heading Chord approximation for exact trajectory mapping
        sh = np.sin(dh * 0.5)
        ch = np.cos(dh * 0.5)
        
        if abs(dh) <= 1e-6:
            sinc = 1.0
        else:
            sinc = (2.0 * sh) / dh
            
        deltaY = dist * sinc
        
        cosH = np.cos(theta_sim[i-1] + dh * 0.5)
        sinH = np.sin(theta_sim[i-1] + dh * 0.5)
        
        # Odom state estimate evolution
        x_sim[i] = x_sim[i-1] + deltaX * cosH + deltaY * sinH
        y_sim[i] = y_sim[i-1] + deltaY * cosH - deltaX * sinH
        theta_sim[i] = anglewrap(theta_sim[i-1] + dh)
        
        lv = lv_new
        
    return x_sim, y_sim

def cost_function(fk, data):
    x_sim, y_sim = simulate(fk, data)
    
    # Target objective: Root Mean Square Error matching vs Particle Filter's final positions
    err = np.sqrt((x_sim - data['pf_x'].values)**2 + (y_sim - data['pf_y'].values)**2)
    return np.mean(err)

def main():
    parser = argparse.ArgumentParser(description="Optimize kinetic friction from driving data.")
    parser.add_argument("--csv", type=str, default="friction_tune.csv", help="Path to the logged CSV file.")
    args = parser.parse_args()

    try:
        data = pd.read_csv(args.csv)
    except FileNotFoundError:
        print(f"File {args.csv} not found! Please make sure you copy it from /usd/friction_tune.csv on the VEX SD Card.")
        return

    if len(data) == 0:
        print("Empty dataset!")
        return

    print(f"Loaded {len(data)} cycles. Optimizing kinetic friction...")
    
    # 1D Nelder-Mead using Brent's method equivalent suitable matching PFConfig values (~68)
    res = minimize_scalar(lambda fk: cost_function(fk, data), bounds=(0, 200), method='bounded')
    
    print("Optimization finished!")
    print(f"Optimal kinetic friction (fk): {res.x:.4f}")
    print(f"Final trajectory loss (mean error): {res.fun:.4f} units")

if __name__ == '__main__':
    main()
