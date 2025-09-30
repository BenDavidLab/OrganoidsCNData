import pandas as pd
import argparse

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter a CNR file to remove rows with log2 < -2.")
    parser.add_argument("cnr_file", help="Path to the input CNR file.")
    parser.add_argument("output_file", help="Path to save the filtered CNR file.")
    return parser.parse_args()

# Main function
def filter_cnr(cnr_file, output_file):
    # Load the CNR file
    cnr_data = pd.read_csv(cnr_file, sep="\t")

    # Filter rows where log2 >= -2
    filtered_cnr = cnr_data[cnr_data["log2"] >= -2]

    # Save the filtered file
    filtered_cnr.to_csv(output_file, sep="\t", index=False)
    print(f"Filtered CNR file saved to: {output_file}")

if __name__ == "__main__":
    args = parse_arguments()
    filter_cnr(args.cnr_file, args.output_file)
