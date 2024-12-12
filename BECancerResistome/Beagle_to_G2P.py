"""
Beagle to G2P Pipeline
September 2024 
Laura Drepanos and Ganna Reint


Purpose: Takes BE tiling data for a single gene and infers codon level effects to visualize in the G2P portal


Example usage: python3 Beagle_to_G2P.py sample_Beagle_to_G2P_input.txt 
Example input file contents:
	Data/MelJuSo_PPM1D_zscores_RDA867.xlsx
	Data/PPM1D-codons.txt
	A-G
	ENST00000305921.8
	sample_output
	
Requirements: Python v.3.8.16, Pandas v.1.5.3


Inputs (provided in input file indicated as command line argument):
1. Path to file containing Z-scores for each sgRNA tiling the gene. 
	- Can be excel, csv, or tab-delimited 
	- Must include following columns: sgRNA Target Sequence, Z-score
2. Path to Per-codon base edits Beagle output file 
3. A-G or C-T depending on whether the screen employed an A-based or C-based editor 
4. Select transcript to use
	- must be present in the "Target Transcript ID" column of codon-level Beagle output
	- ideally canonical MANE Select Transcript
5. Identifier to prepend to output file name 

Outputs: 
1. File (sample_output_Beagle_to_G2P.csv) to upload into G2P portal. Each row is an amino acid along the protein
	Columns: 
	- Amino Acid Position: Integer indicating position of amino acid along the protein
	- Ref AA : Reference Amino Acid at this position on MANE Select Transcript
	*note: next 5 columns are blank if there are no sgRNAs that make their most severe edit at this position,
		since these metrics only consider sgRNAs that generate their most severe edit at this position * 
	- Strongest Z-Score Alt AA: Resulting Amino Acid from Strongest Z-Score sgRNA
	- Strongest Z-Score sgRNA: sgRNA with largest (absolute value) z-score that targets this codon
	- Strongest Z-Score: Z-score of ^
	- Strongest Z-Score mutation: Mutation introduced by ^^ 
	- Mean Z-Score: Mean Z-score of all sgRNAs targeting this codon 
	- Most Severe Possible Mutation: The most severe mutation type that can occur at this position  
		- DOES include sgRNAs that generate a more severe edit at another codon
		- not necessarily the sgRNA with the strongest Z-score targeting this codon 

"""

import sys
import pandas as pd


def load_input():
    if len(sys.argv) < 2:
        print(
            "\nPlease provide filepath for inputs\nExample usage:\npython3 Beagle_to_G2P.py sample_Beagle_to_G2P_input.txt \n"
        )
        exit()

    # Get input info from input file
    inputfile = open(sys.argv[1], "r")
    inputs = inputfile.readlines()
    inputs = [
        item.replace("\n", "") for item in inputs
    ]  # gets rid of \n character in inputs
    inputfile.close()

    if len(inputs) < 5:
        print(
            "\nPlease ensure input file indicates Z-score file, Beagle output, indication of A-G or C-T editing, desired transcript, and output file identifier\n"
        )
        exit()

    # Read in file with Z-scores for each sgRNA tiling the gene
    zscores_filepath = inputs[0]
    if zscores_filepath[-4:] == "xlsx":
        # if the z-scores file is an excel spreadsheet, read in using read_excel
        # note that this function will, by default, use data from the FIRST sheet of the excel spreadsheet
        zscores = pd.read_excel(zscores_filepath)
    elif zscores_filepath[-3:] == "csv":
        zscores = pd.read_csv(zscores_filepath)
    else:
        zscores = pd.read_table(zscores_filepath)
    # Verify that z-score data includes necessary columns
    if ("sgRNA Target Sequence" not in zscores.columns.tolist()) | (
        "Z-score" not in zscores.columns.tolist()
    ):
        print(
            "\nPlease ensure Z-score file contains the following columns: sgRNA Target Sequence, Z-score\n"
        )
        exit()
    zscores = zscores[["sgRNA Target Sequence", "Z-score"]]

    # Read in Beagle codon-level output
    beagle_filepath = inputs[1]
    if beagle_filepath[-3:] != "txt":
        print(
            "\nPlease ensure second input is the path to Per Codon base edits .txt file downloaded from Beagle"
        )
        exit()
    beagle_data = pd.read_table(beagle_filepath)

    # Identify if the tiling screen employed A-based editing or C-based editing
    edit_type = inputs[2]
    if (edit_type != "A-G") and (edit_type != "C-T"):
        print(
            "\nThird line of input file should be A-G or C-T depending on the base editor employed\n"
        )
        exit()

    # Indicate transcript to use
    transcript = inputs[3]
    if transcript not in beagle_data["Target Transcript ID"].tolist():
        print(
            "\nFourth line of input file must indicate a transcript that is present in 'Target Transcript ID' column of Beagle codon-level output \n"
        )
        exit()

    # Get identifier for output file name
    identifier = inputs[4]

    return zscores, beagle_data, edit_type, transcript, identifier


# Get sgRNA-level base-editing data relative to amino acid
def guide_to_residue(zscores, beagle_data, edit_type, transcript):

    # Obtain reference and mutant amino acid from variant notation
    beagle_data["Ref AA"] = beagle_data["Amino Acid Edit"].str[:3]
    beagle_data["Alt AA"] = beagle_data["Amino Acid Edit"].str[-3:]

    # From Beagle output, filter out edits that do not pertain to editor or target transcript or are empty
    beagle_data = beagle_data[
        (beagle_data["Edit Type"] == edit_type)
        & (beagle_data["Target Transcript ID"] == transcript)
    ].reset_index(drop=True)
    beagle_data = beagle_data[-beagle_data["Mutation Category"].isna()].reset_index(
        drop=True
    )

    # Merge Z-score and Beagle output
    merged_data = beagle_data.merge(zscores, on="sgRNA Target Sequence")

    # Limit bystander edits by only retaining most severe edit predicted by each guide

    severity_rank = {
        "Splice-donor": 1,
        "Splice-acceptor": 2,
        "Nonsense": 3,
        "Missense": 4,
        "Silent": 5,
        "UTR": 6,
        "Intron": 7,
        "Flank": 8,
        "Unknown": 9,
    }
    merged_data["Mutation Severity"] = merged_data["Mutation Category"].apply(
        lambda x: severity_rank[x]
    )
    # Identify highest severity mutation caused by each guide
    sgRNA_top_severity_df = (
        merged_data.groupby(["sgRNA Target Sequence"])
        .agg(sgRNA_top_severity=("Mutation Severity", "min"))
        .reset_index()
    )
    merged_data_scored = merged_data.merge(
        sgRNA_top_severity_df, on="sgRNA Target Sequence"
    )
    # For each codon, remove edits from guides that make a more severe edit at another codon
    merged_data_topseverity = (
        merged_data_scored[
            merged_data_scored["Mutation Severity"]
            == merged_data_scored["sgRNA_top_severity"]
        ]
        .reset_index(drop=True)
        .copy()
    )

    # Represent each codon with one row only
    condense_by_pos = (
        merged_data_topseverity.groupby(["Amino Acid Position", "Ref AA"])
        .agg(
            max_z=("Z-score", "max"),
            min_z=("Z-score", "min"),
            mean_z=("Z-score", "mean"),
        )
        .reset_index()
    )

    # Retain the z-score with the highest absolute value at each codon
    condense_by_pos["strongest_z"] = condense_by_pos[["max_z", "min_z"]].apply(
        lambda x: max(x.min(), x.max(), key=abs), axis=1
    )

    # Rejoin condensed data with full dataset to associate the strongest z-score edit with the corresponding sgRNA sequence, amino acid identity, and mutation type
    full_with_best_z = condense_by_pos.merge(
        merged_data_topseverity, on=["Amino Acid Position", "Ref AA"]
    )
    full_with_best_z = (
        full_with_best_z[full_with_best_z["Z-score"] == full_with_best_z["strongest_z"]]
        .drop_duplicates(keep="first")
        .reset_index(drop=True)
    )
    codon_level_with_best_z = full_with_best_z[
        [
            "Amino Acid Position",
            "Ref AA",
            "mean_z",
            "strongest_z",
            "Alt AA",
            "sgRNA Target Sequence",
            "Mutation Category",
            "Mutation Severity",
        ]
    ]

    # Identify the most severe mutation level possible at each codon. May be due to sgRNA that makes a more severe edit at another codon
    highest_severity_by_residue = (
        merged_data_scored.groupby(["Amino Acid Position", "Ref AA"])
        .agg(max_severity=("Mutation Severity", "min"))
        .reset_index()
    )
    max_severity_with_max_z_score = highest_severity_by_residue.merge(
        codon_level_with_best_z, on=["Amino Acid Position", "Ref AA"], how="left"
    )
    # Associate mutation severity level (numeric) back to mutation type
    max_severity_with_max_z_score["Most Severe Possible Mutation"] = (
        max_severity_with_max_z_score["max_severity"].apply(
            lambda x: list(severity_rank.keys())[list(severity_rank.values()).index(x)]
        )
    )

    return max_severity_with_max_z_score


# Gets processed data into file compatible with G2P
def to_output(processed_data, outputfilename):

    # subset and rename columns to be informative output
    g2p_input = processed_data[
        [
            "Amino Acid Position",
            "Ref AA",
            "Alt AA",
            "sgRNA Target Sequence",
            "strongest_z",
            "Mutation Category",
            "mean_z",
            "Most Severe Possible Mutation",
        ]
    ]
    g2p_input = g2p_input.rename(
        columns={
            "Alt AA": "Strongest Z-Score Alt AA",
            "sgRNA Target Sequence": "Strongest Z-Score sgRNA",
            "strongest_z": "Strongest Z-score",
            "Mutation Category": "Strongest Z-Score mutation",
            "mean_z": "Mean Z-Score",
        }
    )

    output_filename = outputfilename + "_Beagle_to_G2P.csv"
    g2p_input.to_csv(output_filename, index=False)
    print("file created:", output_filename)


def main():
    zscores, beagle_data, edit_type, transcript, outputfilename = load_input()
    processed_data = guide_to_residue(zscores, beagle_data, edit_type, transcript)
    to_output(processed_data, outputfilename)


if __name__ == "__main__":
    main()
