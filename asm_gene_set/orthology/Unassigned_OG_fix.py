import pandas as pd
import sys

def main(input_file):
    # Read TSV file
    df = pd.read_csv(input_file, sep='\t', dtype=str).fillna("")

    # Melt to long format: each row will be Orthogroup, Source, Gene
    df_long = df.melt(id_vars=[df.columns[0]], var_name="Source", value_name="Gene")

    # Filter out empty entries
    df_filtered = df_long[df_long["Gene"] != ""]

    # Output the result as TSV
    output_file = input_file.replace(".tsv", "_gene_assignments.tsv")
    df_filtered.to_csv(output_file, sep='\t', index=False)
    print(f"Output written to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python extract_gene_assignments.py <input_file.tsv>")
        sys.exit(1)

    input_file = sys.argv[1]
    main(input_file)
