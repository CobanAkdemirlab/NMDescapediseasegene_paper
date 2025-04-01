import metapredict as meta
import os

# Directory to save the results
output_dir = "plus2_disorder_scores"  # Directory to save the results

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Loop through all FASTA files in the current directory
for fasta_file in os.listdir():  # List files in the current directory
    if fasta_file.endswith(".fasta") or fasta_file.endswith(".fa"):  # Check for FASTA files
        fasta_path = os.path.join(os.getcwd(), fasta_file)  # Get absolute path

        # Read the sequence from the FASTA file
        with open(fasta_path, "r") as file:
            lines = file.readlines()
            test_seq = "".join(line.strip() for line in lines if not line.startswith(">"))  # Remove header

        # Predict disorder scores
        test_scores = meta.predict_disorder(test_seq)

        # Save the results to a file
        output_file = os.path.join(output_dir, f"{os.path.splitext(fasta_file)[0]}_disorder_scores.txt")
        with open(output_file, "w") as file:
            file.write(f"Sequence: {test_seq}\n")
            file.write("Disorder Scores:\n")
            file.write(", ".join(map(str, test_scores)))  # Convert scores to a comma-separated string

        print(f"Results saved to {output_file}")

