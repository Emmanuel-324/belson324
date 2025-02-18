import numpy as np

# Define the given plane normal and direction
plane_normal = np.array([1, 1, 0], dtype=float)  # (312) plane normal
direction = np.array([0, 0, 1], dtype=float)  # [021̅] crystal direction

# Normalize the given direction to get the new x-axis
x_new = direction / np.linalg.norm(direction)

# Compute the new y-axis using the cross product (to ensure orthogonality)
y_new = np.cross(plane_normal, direction).astype(float)  
y_new /= np.linalg.norm(y_new)  # Normalize

# The new z-axis is the plane normal, normalized
z_new = plane_normal / np.linalg.norm(plane_normal)

# Construct the rotation matrix
R = np.column_stack((x_new, y_new, z_new))

# Print the rotation matrix
print("Rotation matrix for (312)[021̅]:")
print(R)



import numpy as np

# Define Euler angles in degrees
euler_angles = [0, 0, 45]  # (phi1, Phi, phi2) in degrees

# Convert degrees to radians
euler_angles = np.radians(euler_angles)
phi1, Phi, phi2 = euler_angles

# Compute the corrected Bunge rotation matrix
R = np.array([
    [np.cos(phi1) * np.cos(phi2) - np.sin(phi1) * np.sin(phi2) * np.cos(Phi),
     np.cos(phi1) * np.sin(phi2) + np.sin(phi1) * np.cos(phi2) * np.cos(Phi),
     np.sin(phi2) * np.sin(Phi)],

    [-np.sin(phi1) * np.cos(phi2) - np.cos(phi1) * np.sin(phi2) * np.cos(Phi),
     -np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(Phi),
     np.cos(phi2) * np.sin(Phi)],

    [np.sin(phi1) * np.sin(Phi),
     -np.cos(phi1) * np.sin(Phi),
     np.cos(Phi)]
])

# Standard basis vectors
x_std = np.array([1, 0, 0])  # Standard x-direction
y_std = np.array([0, 1, 0])  # Standard y-direction
z_std = np.array([0, 0, 1])  # Standard z-direction

# Compute the transformed crystal directions
hkl = R @ z_std  # New z-axis (Miller indices of the plane normal)
uvw = R @ x_std  # New x-axis (Miller indices of the crystal direction)

# Round to nearest integers with a correction for small floating-point errors
hkl = np.round(hkl).astype(int)
uvw = np.round(uvw).astype(int)

# Ensure correct Miller notation for values close to -1
uvw[np.isclose(uvw, -1, atol=1e-3)] = -1
uvw[np.isclose(uvw, 1, atol=1e-3)] = 1
uvw[np.isclose(uvw, 0, atol=1e-6)] = 0  # Ensure very small values become zero

# Print the results
print("Miller indices for Euler angles (0, 45,0):")
print(f"(hkl) = ({hkl[0]}, {hkl[1]}, {hkl[2]})")
print(f"[uvw] = [{uvw[0]}, {uvw[1]}, {uvw[2]}]")
