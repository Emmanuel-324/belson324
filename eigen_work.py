import pandas as pd

def process_csv(file_path, column1, column2, output_file):
    """
    Subtracts column2 from column1 in a CSV file, takes the absolute value,
    and saves the modified data back as a CSV file.

    Parameters:
    file_path (str): Path to the input CSV file.
    column1 (str): Name of the first column.
    column2 (str): Name of the second column.
    output_file (str): Path to save the modified CSV file.
    """
    # Load the CSV file
    df = pd.read_csv(file_path)

    # Ensure the specified columns exist
    if column1 not in df.columns or column2 not in df.columns:
        raise ValueError(f"Columns '{column1}' or '{column2}' not found in the file.")

    # Subtract and take absolute value
    df["Result"] = abs(df[column1] - df[column2])

    # Save to a new CSV file
    df.to_csv(output_file, index=False)
    print(f"Processed file saved as: {output_file}")

# Example usage
file_path = "/home/emmanuel324/projects/belson324/Crystal_Plasticity/2-14-25/Increase_eigen/Test1/Data_eigen_increase/eigen_increase_out.csv"  # Your actual file path
column1 = "eth_xx"  # Replace with your actual column name
column2 = "strain_xx"  # Replace with your actual column name
output_file = file_path  # Save in the same file

process_csv(file_path, column1, column2, output_file)

