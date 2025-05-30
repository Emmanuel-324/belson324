

import pandas as pd

# Load the CSV file
input_file = "/home/emmanuel324/projects/belson324/Crystal_Plasticity/2-28-25/Eigen/parameters1_eigen_new/Data/parameters1_eigen_new_out.csv"
data = pd.read_csv(input_file)
df = pd.read_csv(input_file)
# Ensure 'time' column 
if 'time' not in data.columns:
    raise ValueError("The input CSV file must contain a 'time' column")


L = 300  #length of domain

df['stress_xx'] = df['stress_xx'].abs()

# Calculate displacement (delta) and strain
data['displacement'] = data['time'] * 0.1    
data['strain'] = data['displacement'] / L

# Save the updated data to a new CSV file
output_file = input_file
data.to_csv(output_file, index=False)

print(f"Output with strain saved to {output_file}")


