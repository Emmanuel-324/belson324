# Read the input slip system file
with open('input_slip_sys (4).txt', 'r') as file:
    slip_systems = file.readlines()

# Define the selection criteria for slip planes and directions
desired_planes = [(1, 1, 1), (1, -1, 1), (-1, 1, 1)]  # Example: {111} planes
desired_directions = [(1, 1, 0), (1, -1, 0)]  # Example: <110> directions

# Filter slip systems based on the criteria
selected_slip_systems = []

for system in slip_systems:
    plane_normal = tuple(map(int, system.split()[:3]))  # First three values: plane normal
    direction = tuple(map(int, system.split()[3:]))  # Last three values: direction

    if plane_normal in desired_planes and direction in desired_directions:
        selected_slip_systems.append(system)

# Print the selected slip systems
for system in selected_slip_systems:
    print(system)
