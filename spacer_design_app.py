import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from io import StringIO, BytesIO

try:
    import streamlit as st
except ModuleNotFoundError:
    raise ModuleNotFoundError("The 'streamlit' module is not installed. Please install it using 'pip install streamlit'.")

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

# Upload the GenBank file with tooltip
genbank_file = st.file_uploader("Upload a GenBank file (.gb or .gbff) ❗", type=["gb", "gbff"], help="Upload the GenBank file containing your target gene's sequence.")

# User inputs for parameters
target_genes = st.text_input("Target Genes (comma-separated) ❗", value="edd", help="Enter the gene names for which you want to design spacers. Use commas to separate multiple genes.")
pam_sequence = st.text_input("PAM Sequence ❗", value="CC", help="Specify the PAM sequence for your chosen CRISPR system (e.g., NGG for SpCas9).")
spacer_length = st.slider("Spacer Length ❗", min_value=16, max_value=50, value=32, help="Select the length of the spacer sequence.")
direction = st.selectbox("Spacer Direction ❗", ["PAM-Upstream (PAM-Spacer)", "PAM-Downstream (Spacer-PAM)"], help="Choose whether the PAM appears upstream or downstream of the spacer.")
design_upstream = st.checkbox("Design Spacers Upstream of ORF ❗", help="Check this box if you want to design spacers upstream of the ORF instead of within the gene.")
upstream_window = st.slider("Upstream Window for PAM (bases from ORF start) ❗", min_value=10, max_value=500, value=(10, 100), help="Select the range (in bases) upstream of the ORF where the PAM should be located.")
num_spacers = st.slider("Number of Spacers to Predict per Gene ❗", min_value=1, max_value=10, value=3, help="Choose how many spacers to predict per gene.")
position_range = st.slider("Search Position in Gene (% of Gene Length) ❗", min_value=0, max_value=100, value=(10, 20), help="Select the percentage range within the gene to search for spacers.")
left_flank = st.text_input("Left Flank Sequence ❗", value="aggtcTcaaaac", help="Specify the left flank sequence to be added before the spacer.")
right_flank = st.text_input("Right Flank Sequence ❗", value="gtttttGAGACCa", help="Specify the right flank sequence to be added after the spacer.")

# Function to extract ORF start sites and upstream sequences from GenBank
def get_orf_upstream_sequences(gb_content, gene_names, upstream_window):
    gene_sequences = {}
    upstream_sequences = {}
    gb_file = StringIO(gb_content)
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
                if gene_name in gene_names:
                    orf_start = feature.location.start
                    strand = feature.location.strand
                    gene_seq = feature.location.extract(record).seq
                    gene_sequences[gene_name] = gene_seq
                    
                    # Extract upstream region considering strand orientation
                    if strand == 1:  # Forward strand
                        upstream_start = max(0, orf_start - upstream_window[1])
                        upstream_end = max(0, orf_start - upstream_window[0])
                        upstream_seq = record.seq[upstream_start:upstream_end]
                    else:  # Reverse strand
                        upstream_start = min(len(record.seq), feature.location.end + upstream_window[0])
                        upstream_end = min(len(record.seq), feature.location.end + upstream_window[1])
                        upstream_seq = record.seq[upstream_start:upstream_end].reverse_complement()
                    
                    upstream_sequences[gene_name] = upstream_seq
    return gene_sequences, upstream_sequences

# Run analysis when the button is clicked
if st.button("Generate Spacers ❗", help="Click to generate spacers based on the selected parameters."):
    if genbank_file is None:
        st.error("Please upload a GenBank file.")
    else:
        try:
            # Read and parse the GenBank file
            gb_content = genbank_file.getvalue().decode("utf-8")
            target_genes_list = [gene.strip() for gene in target_genes.split(",")]

            # Extract ORF and upstream sequences
            gene_sequences, upstream_sequences = get_orf_upstream_sequences(gb_content, target_genes_list, upstream_window)

            # Check if any gene sequences were found
            if not gene_sequences:
                st.error("No target genes found in the GenBank file. Please check the gene names.")
            else:
                # Generate spacers
                spacers = generate_spacers_with_pam(gene_sequences, upstream_sequences, pam_sequence, spacer_length, direction, num_spacers, position_range, design_upstream)

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
                    label="Download Results as CSV ❗", data=output, file_name="spacers_output.csv", mime="text/csv", help="Download the generated spacer sequences as a CSV file.")

        except Exception as e:
            st.error(f"An error occurred: {e}")
