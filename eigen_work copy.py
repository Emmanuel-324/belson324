import pandas as pd

def process_stress_strain(csv_file):
  
    # Load the data
    df = pd.read_csv(csv_file)

    # Ensure necessary columns exist
    if "strain_xx" not in df.columns or "stress_xx" not in df.columns:
        raise ValueError("CSV file must contain 'strain_xx' and 'stress_xx' columns.")

    # Extract the tensile region (stress > 0)
    df_tensile = df[df["stress_xx"] > 0].copy()

    # Shift the strain values so the first tensile strain starts at zero
    df_tensile["strain_xx"] = df_tensile["strain_xx"] - df_tensile["strain_xx"].iloc[0]

    # Save the processed data back to the same CSV file
    df_tensile.to_csv(csv_file, index=False)
    print(f"Processed tensile stress-strain data saved to: {csv_file}")

# Example usage:
csv_file_path = "/home/emmanuel324/projects/belson324/Crystal_Plasticity/2-14-25/Slip_Modes/slip_mode2/Data_slip_system_mode2/eigen_increase_out.csv"  # Replace this with the actual file path
process_stress_strain(csv_file_path)
