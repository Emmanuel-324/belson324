import pandas as pd

# Load the CSV file with error handling for malformed rows
input_file = "/home/emmanuel324/projects/belson/Crystal_Plasticity_Test/Data/Test_s184_out.csv"
try:
    # Attempt to read the CSV file and skip bad lines
    data = pd.read_csv(input_file, on_bad_lines='skip')
    print(f"CSV file loaded successfully. Number of rows: {len(data)}")
except pd.errors.ParserError as e:
    print(f"Error reading CSV file: {e}")
    raise

# Ensure 'time' column exists
if 'time' not in data.columns:
    raise ValueError("The input CSV file must contain a 'time' column")

# Convert 'time' column to numeric, coercing errors to NaN for non-numeric values
data['time'] = pd.to_numeric(data['time'], errors='coerce')

# Drop rows with NaN values in 'time' column (if any)
data = data.dropna(subset=['time'])

# Verify that 'time' is now numeric and print sample values
print(f"'time' column data type: {data['time'].dtype}")
print(f"First few entries in 'time' column:\n{data['time'].head()}")

# Define the domain length
L = 300  # Replace with the actual length of your domain

# Calculate displacement and strain
try:
    # Ensure 'time' is numeric before performing arithmetic operations
    data['displacement'] = data['time'].astype(float) * 0.1
    data['strain'] = data['displacement'] / L
except Exception as e:
    print(f"Error occurred during calculation: {e}")
    raise

# Save the updated data to a new CSV file
output_file = "/home/emmanuel324/projects/belson/Crystal_Plasticity_Test/Data/Test_s184_out.csv"
data.to_csv(output_file, index=False)

print(f"Output with strain saved to {output_file}")
