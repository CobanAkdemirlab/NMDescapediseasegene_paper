import os
import glob
import metapredict as meta
import pandas as pd

# Define the directory containing FASTA files
fasta_files = glob.glob("*.fasta")

# Output storage
results = []

# Function to read sequences from FASTA files
def read_fasta(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        header = lines[0].strip()[1:]  # Remove ">" from the header
        sequence = "".join(line.strip() for line in lines[1:])
    return header, sequence

# Process each FASTA file
for fasta_file in fasta_files:
    key, sequence = read_fasta(fasta_file)

    # Predict disorder
    disorder_domains = meta.predict_disorder_domains(sequence)

    # Store results
    results.append({
        "Variant_Key": key,
        "Disorder_Domains": disorder_domains.disordered_domain_boundaries
    })

    # Print summary
    print(f"Processed {fasta_file}:")
    print(f"  Disorder Domains: {disorder_domains}")

# Save results to a CSV file
df = pd.DataFrame(results)
df.to_csv("disorder_predictions.csv", index=False)

print("Disorder predictions saved to disorder_predictions.csv")
