import numpy as np
import matplotlib.pyplot as plt

class LidarSimulator:
    def __init__(self, divergence_deg=2.0, num_rays=11, aperture_w0=0.01):
        """
        Simulates a LiDAR sensor with beam divergence.
        :param divergence_deg: Total divergence angle in degrees.
        :param num_rays: Number of rays in the bundle.
        :param aperture_w0: Initial width of the laser aperture.
        """
        self.theta_div = np.radians(divergence_deg)
        self.num_rays = num_rays
        self.w0 = aperture_w0
        
    def get_ray_bundle(self, heading_deg):
        """Returns a list of angles for the rays in the bundle."""
        center_angle = np.radians(heading_deg)
        half_div = self.theta_div / 2.0
        return np.linspace(center_angle - half_div, center_angle + half_div, self.num_rays)

    def calculate_gaussian_weight(self, ray_angle, center_angle):
        """
        Calculates the Gaussian weight for a ray.
        Weight_i = e^(-2 * phi_i^2 / (theta_div / 2)^2)
        """
        phi = ray_angle - center_angle
        # Normalizing by half-divergence
        sigma = self.theta_div / 2.0
        if sigma == 0: return 1.0
        return np.exp(-2 * (phi**2) / (sigma**2))

    def cast_ray(self, origin, angle, walls):
        """
        Finds the nearest intersection of a ray with defined walls.
        Returns (dist, hit_point, hit_normal)
        """
        min_dist = float('inf')
        hit_point = None
        hit_normal = None
        
        ox, oy = origin
        dx, dy = np.cos(angle), np.sin(angle)
        
        for (p1, p2) in walls:
            x1, y1 = p1
            x2, y2 = p2
            vx, vy = x2 - x1, y2 - y1
            
            det = vy * dx - vx * dy
            if abs(det) < 1e-9:
                continue
                
            t = (vx * (oy - y1) - vy * (ox - x1)) / det
            u = (dx * (oy - y1) - dy * (ox - x1)) / det
            
            if t > 0 and 0 <= u <= 1:
                if t < min_dist:
                    min_dist = t
                    hit_point = (ox + t * dx, oy + t * dy)
                    # Normal is perpendicular to (vx, vy). Point towards sensor.
                    n = np.array([-vy, vx], dtype=float)
                    mag = np.sqrt(vx**2 + vy**2)
                    if mag > 0: n = n / mag
                    if np.dot(n, [dx, dy]) > 0: n = -n
                    hit_normal = n
                    
        return min_dist, hit_point, hit_normal

    def simulate_reading(self, origin, heading_deg, walls):
        """
        Simulates the sensor reading using a bundle of rays.
        Returns (final_distance_strongest, final_distance_blurred, best_ray_angle, ray_data)
        """
        ray_angles = self.get_ray_bundle(heading_deg)
        center_angle = np.radians(heading_deg)
        
        results = []
        for angle in ray_angles:
            dist, hit, normal = self.cast_ray(origin, angle, walls)
            weight = self.calculate_gaussian_weight(angle, center_angle)
            
            if normal is not None:
                ray_dir = np.array([np.cos(angle), np.sin(angle)])
                cos_alpha = abs(np.dot(normal, -ray_dir))
            else:
                cos_alpha = 0.0
            
            power = (weight * cos_alpha) / (dist**2) if dist < float('inf') and dist > 0.01 else 0
            
            results.append({
                'angle': angle,
                'dist': dist,
                'hit': hit,
                'normal': normal,
                'weight': weight,
                'cos_alpha': cos_alpha,
                'power': power,
            })
            
        powers = [r['power'] for r in results]
        if max(powers) > 0:
            strongest_idx = np.argmax(powers)
            d_strongest = results[strongest_idx]['dist']
            best_ray_angle = results[strongest_idx]['angle']
        else:
            d_strongest = float('inf')
            best_ray_angle = center_angle
        
        total_power = sum(powers)
        if total_power > 0:
            valid_results = [r for r in results if r['dist'] < float('inf')]
            d_blurred = sum(r['dist'] * r['power'] for r in valid_results) / total_power
        else:
            d_blurred = float('inf')
            
        print(f"Strongest: {d_strongest}m (from ray {np.degrees(best_ray_angle):.1f} deg)")
            
        return d_strongest, d_blurred, best_ray_angle, results

def plot_simulation(origin, heading, walls, simulator):
    d_strong, d_blur, best_angle, ray_data = simulator.simulate_reading(origin, heading, walls)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    ax1.set_title(f"Distance Sensor Beam Divergence (Heading: {heading} deg)")
    for i, (p1, p2) in enumerate(walls):
        ax1.plot([p1[0], p2[0]], [p1[1], p2[1]], 'k-', lw=2, label="Wall" if i==0 else "")
    
    ax1.plot(origin[0], origin[1], 'ro', label="Sensor")
    
    for r in ray_data:
        if r['hit']:
            ax1.plot([origin[0], r['hit'][0]], [origin[1], r['hit'][1]], 
                     color='orange', alpha=0.3, lw=1)
            ax1.scatter(r['hit'][0], r['hit'][1], color='orange', s=20 * r['weight'], alpha=0.6)
            
    h_rad = np.radians(heading)
    
    ax1.plot([origin[0], origin[0] + 5 * np.cos(h_rad)],
             [origin[1], origin[1] + 5 * np.sin(h_rad)], 
             'r-.', lw=1, alpha=0.5, label=f"Heading ({heading}deg)")

    if d_strong < float('inf'):
        ax1.plot([origin[0], origin[0] + d_strong * np.cos(best_angle)],
                 [origin[1], origin[1] + d_strong * np.sin(best_angle)], 
                 'm-', lw=1, alpha=0.8, label=f"Raw Strongest Ray ({np.degrees(best_angle):.1f}deg)")

    if d_strong < float('inf'):
        ax1.plot([origin[0], origin[0] + d_strong * np.cos(h_rad)],
                 [origin[1], origin[1] + d_strong * np.sin(h_rad)], 
                 'b--', lw=3, label=f"Reported Strongest ({d_strong:.2f}m)")
        
    if d_blur < float('inf'):
        ax1.plot([origin[0], origin[0] + d_blur * np.cos(h_rad)],
                 [origin[1], origin[1] + d_blur * np.sin(h_rad)], 
                 'g:', lw=3, label=f"Blurred ({d_blur:.2f}m)")

    ax1.set_xlim(-1, 5)
    ax1.set_ylim(-1, 5)
    ax1.set_aspect('equal')
    ax1.legend(fontsize='small')
    ax1.grid(True)
    
    ax2.set_title("Signal Power vs Ray Angle")
    angles_relative = [np.degrees(r['angle'] - np.radians(heading)) for r in ray_data]
    powers = [r['power'] for r in ray_data]
    
    ax2.bar(angles_relative, powers, width=0.2, color='orange', alpha=0.7, label="Signal Power (incl. physics)")
    ax2.set_xlabel("Angle Relative to Center (deg)")
    ax2.set_ylabel("Signal Power (Return Signal)")
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig('distance_sensor_sim.png')
    print(f"Simulation saved to 'distance sensor sim.png'")
    plt.show()


if __name__ == "__main__":
    walls = [
        ((0, 0), (4, 0)),
        ((4, 0), (4, 4)),
        ((4, 4), (0, 4)),
        ((0, 4), (0, 0)),
        ((2-.11/2, 4-.12), (2+.11/2, 4-.12))
    ]

    #walls = [((-70.5, -70.5), (70.5, -70.5)),
    #         ((70.5, -70.5), (70.5, 70.5)),
    #         ((70.5, 70.5), (-70.5, 70.5)),
    #         ((-70.5, 70.5), (-70.5, -70.5))]
    
    origin = (2, 4-.510)
    #origin = (39.825, 39.9232)
    #heading = 90 - 44.7
    heading = 90
    # 560
    sim = LidarSimulator(divergence_deg=36.0, num_rays=8)
    
    print(f"Simulating beam hitting a corner from {origin} at {heading} degrees...")
    plot_simulation(origin, heading, walls, sim)
