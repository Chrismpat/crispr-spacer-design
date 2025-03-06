import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from io import StringIO, BytesIO

# Streamlit App Title
st.title("CRISPR Spacer Design Tool")

# Instructions for GBFF File
st.markdown(
    """### How to Download Genome Files from NCBI:
    1. Go to [NCBI Genome](https://www.ncbi.nlm.nih.gov/genome/)
    2. Search for your organism of interest.
    3. Click on the organism name to access its genome page.
    4. Under "Downloads", look for **"GenBank (Full)"**.
    5. Download the **.gbff** file.
    6. Upload the file using the button below.
    
    **Note:** The GBFF file contains both the genome sequence and gene annotations.
    """
)

# Upload the GenBank file
genbank_file = st.file_uploader("Upload a GenBank file (.gb or .gbff)", type=["gb", "gbff"])

# User inputs for parameters
target_genes = st.text_input("Target Genes (comma-separated):", value="edd")
pam_sequence = st.text_input("PAM Sequence:", value="CC")
spacer_length = st.slider("Spacer Length:", min_value=16, max_value=50, value=32)
direction = st.selectbox("Spacer Direction:", ["PAM-Upstream (PAM-Spacer)", "PAM-Downstream (Spacer-PAM)"])
design_upstream = st.checkbox("Design Spacers Upstream of Gene")
upstream_window = st.slider("Upstream Window for PAM (bases from ORF start):", min_value=10, max_value=500, value=(10, 100))
num_spacers = st.slider("Number of Spacers to Predict per Gene:", min_value=1, max_value=10, value=3)
position_range = st.slider("Search Position in Gene (% of Gene Length):", min_value=0, max_value=100, value=(10, 20))
left_flank = st.text_input("Left Flank Sequence:", value="aggtcTcaaaac")
right_flank = st.text_input("Right Flank Sequence:", value="gtttttGAGACCa")

# Function to extract gene sequences from GenBank
def get_gene_sequences(gb_content, gene_names):
    gene_sequences = {}
    gb_file = StringIO(gb_content)
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "gene" and "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
                if gene_name in gene_names:
                    gene_seq = feature.location.extract(record).seq
                    gene_sequences[gene_name] = gene_seq
    return gene_sequences

# Function to generate spacers with PAM
def generate_spacers_with_pam(gene_sequences, pam_seq, spacer_length, direction, num_spacers, position_range, design_upstream, upstream_window):
    spacers = {}
    for gene, sequence in gene_sequences.items():
        gene_spacers = []
        gene_length = len(sequence)
        start_pos = int((position_range[0] / 100) * gene_length)
        end_pos = int((position_range[1] / 100) * gene_length)

        if design_upstream:
            upstream_start = max(0, -upstream_window[1])
            upstream_end = max(0, -upstream_window[0])
        else:
            upstream_start = start_pos
            upstream_end = end_pos

        for i in range(upstream_start, min(upstream_end, gene_length - spacer_length - len(pam_seq))):
            if sequence[i:i + len(pam_seq)] == pam_seq:
                if direction == "PAM-Upstream (PAM-Spacer)":
                    spacer = sequence[i + len(pam_seq):i + len(pam_seq) + spacer_length]
                else:  # "PAM-Downstream (Spacer-PAM)"
                    spacer = sequence[i - spacer_length:i] if i >= spacer_length else ""
                
                if spacer:
                    gene_spacers.append(str(spacer))
            if len(gene_spacers) == num_spacers:
                break
        spacers[gene] = gene_spacers
    return spacers

# Run analysis when the button is clicked
if st.button("Generate Spacers"):
    if genbank_file is None:
        st.error("Please upload a GenBank file.")
    else:
        try:
            # Read and parse the GenBank file
            gb_content = genbank_file.getvalue().decode("utf-8")
            target_genes_list = [gene.strip() for gene in target_genes.split(",")]

            # Extract gene sequences
            gene_sequences = get_gene_sequences(gb_content, target_genes_list)

            # Check if any gene sequences were found
            if not gene_sequences:
                st.error("No target genes found in the GenBank file. Please check the gene names.")
            else:
                # Generate spacers
                spacers = generate_spacers_with_pam(gene_sequences, pam_sequence, spacer_length, direction, num_spacers, position_range, design_upstream, upstream_window)

                # Prepare data for display and download
                spacer_data = []
                for gene, spacer_list in spacers.items():
                    for i, spacer in enumerate(spacer_list):
                        full_spacer = f"{left_flank}{spacer}{right_flank}"
                        spacer_name = f"{gene}_sp.{i + 1}"
                        spacer_data.append([spacer_name, full_spacer])

                # Display results
                df = pd.DataFrame(spacer_data, columns=["Name", "Sequence"])
                st.write("Generated Spacers:")
                st.dataframe(df)

                # Provide CSV download option
                output = BytesIO()
                df.to_csv(output, index=False)
                output.seek(0)
                st.download_button(
                    label="Download Results as CSV",
                    data=output,
                    file_name="spacers_output.csv",
                    mime="text/csv"
                )

        except Exception as e:
            st.error(f"An error occurred: {e}")
