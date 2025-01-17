import numpy as np

def normalize(v):
    return v / np.linalg.norm(v)

def rotation_matrix_to_euler(R):
    """
    Convert rotation matrix to Euler angles in the Bunge convention.
    """
    if R[2, 2] < 1:
        if R[2, 2] > -1:
            phi1 = np.arctan2(R[1, 2], R[0, 2])
            Phi = np.arccos(R[2, 2])
            phi2 = np.arctan2(R[2, 1], -R[2, 0])
        else:  # R[2,2] == -1
            phi1 = -np.arctan2(-R[1, 0], R[1, 1])
            Phi = np.pi
            phi2 = 0
    else:  # R[2,2] == 1
        phi1 = np.arctan2(-R[1, 0], R[1, 1])
        Phi = 0
        phi2 = 0
    return np.degrees(phi1), np.degrees(Phi), np.degrees(phi2)

# Define crystal frame slip direction and plane normal
slip_direction = normalize(np.array([1, 1, 1]))
plane_normal = normalize(np.array([1, 1, 0]))

# Define the global reference frame
z_global = np.array([0, 0, 1])
x_global = np.array([1, 0, 0])

# Rotation to align slip direction with z-axis
axis = np.cross(slip_direction, z_global)
angle = np.arccos(np.dot(slip_direction, z_global))
R1 = np.eye(3)
if np.linalg.norm(axis) > 1e-6:
    axis = normalize(axis)
    ux, uy, uz = axis
    c, s = np.cos(angle), np.sin(angle)
    R1 = np.array([
        [c + ux**2 * (1 - c), ux * uy * (1 - c) - uz * s, ux * uz * (1 - c) + uy * s],
        [uy * ux * (1 - c) + uz * s, c + uy**2 * (1 - c), uy * uz * (1 - c) - ux * s],
        [uz * ux * (1 - c) - uy * s, uz * uy * (1 - c) + ux * s, c + uz**2 * (1 - c)]
    ])

# Adjust plane normal to the correct orientation
adjusted_normal = R1 @ plane_normal

# Final Euler angles
phi1, Phi, phi2 = rotation_matrix_to_euler(R1)
print(f"Euler angles (φ1, Φ, φ2): {phi1:.2f}°, {Phi:.2f}°, {phi2:.2f}°")
