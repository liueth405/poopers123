import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_wtf(file_path='wtf.txt'):
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found.")
        return

    # Read data
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return

    # Convert time to seconds for plotting
    df['t_s'] = (df['t_ms'] - df['t_ms'].iloc[0]) / 1000.0

    # Create figure
    fig = plt.figure(figsize=(15, 10))
    
    # 1. 2D Spatial Plot (Y vs X)
    ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2, rowspan=2)
    ax1.plot(df['ref_x_m'], df['ref_y_m'], 'r--', label='Reference Path', alpha=0.8)
    ax1.plot(df['mpc_x_m'], df['mpc_y_m'], 'b-', label='MPC Actual', linewidth=2)
    
    # Add orientation arrows occasionally
    skip = max(1, len(df)//20)
    for i in range(0, len(df), skip):
        # Reference arrows
        ax1.arrow(df['ref_x_m'].iloc[i], df['ref_y_m'].iloc[i], 
                 0.05 * np.sin(df['ref_th_rad'].iloc[i]), 
                 0.05 * np.cos(df['ref_th_rad'].iloc[i]), 
                 head_width=0.02, color='red', alpha=0.3)
        # MPC arrows
        ax1.arrow(df['mpc_x_m'].iloc[i], df['mpc_y_m'].iloc[i], 
                 0.05 * np.sin(df['mpc_th_rad'].iloc[i]), 
                 0.05 * np.cos(df['mpc_th_rad'].iloc[i]), 
                 head_width=0.02, color='blue', alpha=0.5)

    ax1.set_xlabel('X (meters)')
    ax1.set_ylabel('Y (meters)')
    ax1.set_title('Spatial Trajectory (Field View)')
    ax1.legend()
    ax1.axis('equal')
    ax1.grid(True, linestyle=':', alpha=0.7)

    # 2. X position vs Time
    ax2 = plt.subplot2grid((2, 3), (0, 2))
    ax2.plot(df['t_s'], df['ref_x_m'], 'r--', label='Ref')
    ax2.plot(df['t_s'], df['mpc_x_m'], 'b-', label='MPC')
    ax2.set_ylabel('X (m)')
    ax2.set_title('X vs Time')
    ax2.grid(True)

    # 3. Y position vs Time
    ax3 = plt.subplot2grid((2, 3), (1, 2))
    ax3.plot(df['t_s'], df['ref_y_m'], 'r--', label='Ref')
    ax3.plot(df['t_s'], df['mpc_y_m'], 'b-', label='MPC')
    ax3.set_ylabel('Y (m)')
    ax3.set_xlabel('Time (s)')
    ax3.set_title('Y vs Time')
    ax3.grid(True)
    
    plt.tight_layout()
    
    # Second Figure: Heading and Diagnostics
    fig2, (ax4, ax5) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    # Heading plot (wrapped handling)
    ax4.plot(df['t_s'], np.degrees(df['ref_th_rad']), 'r--', label='Ref Heading')
    ax4.plot(df['t_s'], np.degrees(df['mpc_th_rad']), 'b-', label='MPC Heading')
    ax4.set_ylabel('Heading (deg)')
    ax4.set_title('Heading Tracking')
    ax4.legend()
    ax4.grid(True)
    
    # Cross track error approximation or Vs
    ax5.plot(df['t_s'], df['u_vs_mps'], 'g-', label='Path Progress Speed (Vs)')
    ax5.plot(df['t_s'], df['u_v_mps'], 'k-', label='Robot Velocity (v)', alpha=0.6)
    ax5.set_ylabel('Velocity (m/s)')
    ax5.set_xlabel('Time (s)')
    ax5.set_title('MPC Control Inputs / Progress')
    ax5.legend()
    ax5.grid(True)

    plt.tight_layout()
    print("Graphs generated. Showing plots...")
    plt.show()

if __name__ == "__main__":
    plot_wtf()
