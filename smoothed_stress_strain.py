import pandas as pd

def smooth_stress_strain(csv_file, stress_col="stress_xx", strain_col="strain_xx", threshold=1e-3):
    # Load the CSV file
    df = pd.read_csv(csv_file)
    
    # Remove negative stress values
    df_filtered = df[df[stress_col] >= 0].reset_index(drop=True)
    
    # Identify the initial region where stress is nearly constant
    def find_transition_index(df):
        for i in range(1, len(df)):
            if abs(df[stress_col].iloc[i] - df[stress_col].iloc[i-1]) > threshold:
                return i
        return 0  # If no significant change, return the first index
    
    # Find the transition index
    transition_idx = find_transition_index(df_filtered)
    
    # Remove initial near-constant stress region
    df_smooth = df_filtered.iloc[transition_idx:].reset_index(drop=True)
    
    # Save the smoothed data back to the same CSV file
    df_smooth.to_csv(csv_file, index=False)
    
    print(f"Smoothed data saved to {csv_file}")  # Added print statement
    
    return df_smooth

# Example usage:
smoothed_df = smooth_stress_strain("/home/emmanuel324/projects/belson324/Crystal_Plasticity/2-14-25/Increase_eigen/Test1/Data_eigen_increase/eigen_increase_out.csv")
