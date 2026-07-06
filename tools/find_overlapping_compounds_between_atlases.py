import pandas as pd
import argparse
import os

def generate_metabolite_matrix(input_files_list):
    """
    Processes metabolite files and creates a presence/absence matrix.
    """
    key_cols = ['label', 'inchi_key', 'adduct', 'polarity', 'mz', 'rt_peak', 'rt_min', 'rt_max']
    all_data = []

    for file_path in input_files_list:
        if not os.path.exists(file_path):
            print(f"Warning: File {file_path} not found. Skipping.")
            continue

        # Detect delimiter based on extension
        sep = '\t' if file_path.endswith('.tsv') else ','
        df = pd.read_csv(file_path, sep=sep)

        # Handle label/compound_name logic
        if 'label' not in df.columns and 'compound_name' in df.columns:
            df['label'] = df['compound_name']
        
        # Ensure all required columns exist to avoid KeyError
        for col in key_cols:
            if col not in df.columns:
                df[col] = None

        # Subset to required columns and remove duplicates within the same file
        subset = df[key_cols].drop_duplicates()
        
        # Add the filename as a source identifier
        subset['filename'] = os.path.basename(file_path)
        all_data.append(subset)

    if not all_data:
        return pd.DataFrame()

    # Combine all files into one long-form dataframe
    master_df = pd.concat(all_data, ignore_index=True)

    # Pivot the data: 
    # Index = Unique compounds (defined by key_cols)
    # Columns = Filenames
    # Values = 1 if present (count), 0 otherwise
    matrix = master_df.pivot_table(
        index=key_cols, 
        columns='filename', 
        aggfunc='size', 
        fill_value=0
    )

    # Ensure values are binary (1 if present, 0 if not)
    matrix = (matrix > 0).astype(int)

    # Reset index to work with the key columns as regular columns
    matrix = matrix.reset_index()

    # Count how many times each inchi_key + adduct combination appears in the matrix
    inchi_adduct_counts = (
        matrix.groupby(['inchi_key', 'adduct'])
        .size()
        .reset_index(name='inchi_key_and_adduct_count')
    )
    matrix = matrix.merge(inchi_adduct_counts, on=['inchi_key', 'adduct'], how='left')

    # Sort ascending by rt_peak, then descending by inchi_key_and_adduct_count
    matrix = matrix.sort_values(
        by=['rt_peak', 'inchi_key_and_adduct_count'],
        ascending=[True, False]
    )

    return matrix

def main():
    parser = argparse.ArgumentParser(description="Generate metabolite presence/absence matrix.")
    parser.add_argument('--input-files', required=True, help="CSV file containing a list of metabolite file paths.")
    args = parser.parse_args()

    # Read the list of files from the provided CSV
    try:
        with open(args.input_files, 'r') as f:
            # Assuming one file path per line, potentially comma separated
            files = [line.strip().strip('"') for line in f if line.strip()]
            # If the input-files csv is a single-column CSV, this handles it
            if len(files) == 1 and ',' in files[0]:
                files = [x.strip().strip('"') for x in files[0].split(',')]
    except Exception as e:
        print(f"Error reading input files list: {e}")
        return

    result_matrix = generate_metabolite_matrix(files)

    if not result_matrix.empty:
        result_matrix.to_csv('metabolite_matrix.csv', index=False)
    else:
        print("No data processed.")

if __name__ == "__main__":
    main()