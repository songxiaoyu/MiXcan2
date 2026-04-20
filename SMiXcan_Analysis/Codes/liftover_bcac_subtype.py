import os
import pandas as pd
from pyliftover import LiftOver

# Initialize LiftOver once
lo = LiftOver('hg19', 'hg38')

# LiftOver function
def convert_position(row):
    result = lo.convert_coordinate(row['chr'], row['position'])
    return result[0][1] if result else None

# Directory containing input files
input_dir = "/Users/zhusinan/Downloads/adriana"
output_dir = input_dir

input_file = os.path.join(input_dir, "bcac_2020_subtype_zscores_for_s-mixcan.txt")
output_file = os.path.join(output_dir, "bcac_2020_subtype_zscores_for_s-mixcan_hg38.txt")

print(f"Processing {input_file}...")


# Load data
df = pd.read_csv(input_file, delim_whitespace=True, header=0, low_memory=False, on_bad_lines='skip')

# Standardize chr format
df["position"] = df["var_name"].str.split("_", n=2).str[1].astype(int)
df["chr"] = df["var_name"].str.split("_", n=2).str[0].astype(int)
df['chr'] = df['chr'].apply(lambda x: f'chr{x}')
df['position'] = pd.to_numeric(df['position'], errors='coerce')

# Apply liftover
df['POS_hg38'] = df.apply(convert_position, axis=1)
df = df.dropna(subset=['POS_hg38'])

# Save result
df.to_csv(output_file, index=False, na_rep="")

print(f"✓ Saved {len(df)} lifted positions to {output_file}\n")

