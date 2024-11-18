import ipywidgets as widgets
from IPython.display import display, clear_output
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

# Widget for file path input
file_path_input = widgets.Text(
    value='',
    description='GenBank File Path:',
    placeholder='/Users/christosbatianis/Downloads/GCF_000007565.2_ASM756v2_genomic.gbff',
)

# Input widgets for parameters
gene_input = widgets.Text(
    value='edd',
    description='Target Genes:',
    placeholder='Comma-separated gene names (e.g., edd, pgi)',
)

pam_input = widgets.Text(
    value='CC',
    description='PAM Sequence:',
)

spacer_length_input = widgets.IntSlider(
    value=32,
    min=16,
    max=50,
    step=1,
    description='Spacer Length:',
)

left_flank_input = widgets.Text(
    value='aggtcTcaaaac',
    description='Left Flank:',
)

right_flank_input = widgets.Text(
    value='gtttttGAGACCa',
    description='Right Flank:',
)

# Button to run the analysis
run_button = widgets.Button(
    description='Generate Spacers',
    button_style='success',
)

# Output widget to display results
output = widgets.Output()

# Function to extract gene sequences from the GenBank file
def get_gene_sequences(file_path, gene_names):
    gene_sequences = {}
    try:
        with open(file_path, "r") as gb_file:
            for record in SeqIO.parse(gb_file, "genbank"):
                for feature in record.features:
                    if feature.type == "gene" and "gene" in feature.qualifiers:
                        gene_name = feature.qualifiers["gene"][0]
                        if gene_name in gene_names:
                            gene_seq = feature.location.extract(record).seq
                            gene_sequences[gene_name] = gene_seq
        return gene_sequences
    except Exception as e:
        print(f"Error reading GenBank file: {str(e)}")
        return {}

# Function to generate spacers with PAM
def generate_spacers_with_pam(gene_sequences, pam_seq, spacer_length):
    spacers = {}
    for gene, sequence in gene_sequences.items():
        gene_spacers = []
        for i in range(len(sequence) - spacer_length):
            if sequence[i:i + len(pam_seq)] == pam_seq:
                spacer = sequence[i + len(pam_seq):i + len(pam_seq) + spacer_length]
                gene_spacers.append(str(spacer))
            if len(gene_spacers) == 3:  # Limit to 3 spacers per gene
                break
        spacers[gene] = gene_spacers
    return spacers

# Function to handle button click and process data
def on_button_click(b):
    with output:
        clear_output()

        # Get the file path
        file_path = file_path_input.value.strip()
        if not file_path:
            print("Please enter the GenBank file path.")
            return

        # Get input values
        target_genes = [gene.strip() for gene in gene_input.value.split(',')]
        pam_sequence = pam_input.value
        spacer_length = spacer_length_input.value
        left_flank = left_flank_input.value
        right_flank = right_flank_input.value

        # Extract gene sequences
        gene_sequences = get_gene_sequences(file_path, target_genes)

        # Check if any gene sequences were found
        if not gene_sequences:
            print("No target genes found in the GenBank file. Please check the gene names.")
            return
        
        # Generate spacers
        spacers = generate_spacers_with_pam(gene_sequences, pam_sequence, spacer_length)

        # Prepare data for display and Excel export
        spacer_data = []
        for gene, spacer_list in spacers.items():
            for i, spacer in enumerate(spacer_list):
                # Create sense oligo
                full_spacer = f"{left_flank}{spacer}{right_flank}"
                spacer_name_sense = f"{gene}_sp.{i + 1}_Sense"

                # Create antisense oligo (reverse complement)
                reverse_complement = str(Seq(full_spacer).reverse_complement())
                spacer_name_antisense = f"{gene}_sp.{i + 1}_AntiSense"

                # Append to spacer data
                spacer_data.append([spacer_name_sense, full_spacer])
                spacer_data.append([spacer_name_antisense, reverse_complement])

        # Display results in a DataFrame
        df = pd.DataFrame(spacer_data, columns=["Name", "Sequence"])
        display(df)

        # Save results to Excel
        output_filename = "spacers_output.xlsx"
        df.to_excel(output_filename, index=False)
        print(f"Results saved to {output_filename}")

# Assign the button click event
run_button.on_click(on_button_click)

# Display the widgets
display(file_path_input, gene_input, pam_input, spacer_length_input, left_flank_input, right_flank_input, run_button, output)
