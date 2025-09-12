import pandas as pd

# Load the CSV sheet
file_path = '/home/emmanuel324/projects/belson/KKS_large_largescale/Concentration_New/Test5/Data/Conc5_table.csv'  # Replace with your file path
df = pd.read_csv(file_path)  # Use read_csv for CSV files

# Define the total area of the material (300 nm x 300 nm)
total_material_area = 300 * 300  # 90,000 nm²

# Calculate total precipitate area (since there is only gr1area, total precipitate area is gr1area)
df['total_precipitate_area'] = df['gr1area']

# Calculate the area fraction at each time point
df['area_fraction'] = (df['total_precipitate_area'] / total_material_area) * 100

# Filter out rows where the precipitate area decreases compared to the previous time step
df_filtered = df[df['total_precipitate_area'].diff().fillna(0) >= 0]

# Display the updated dataframe with the new area fraction column, showing only increasing areas
print(df_filtered)

# Optionally save the result back to a CSV file
output_file = file_path
df_filtered.to_csv(output_file, index=False)
