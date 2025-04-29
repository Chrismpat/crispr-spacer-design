import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from io import StringIO, BytesIO

# Attempt to import Streamlit for error clarity
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

# File upload
genbank_file = st.file_uploader(
    "Upload a GenBank file (.gb or .gbff) ?", type=["gb", "gbff"],
    help="Upload the GenBank file containing your target gene's sequence."
)

# User inputs
target_genes = st.text_input(
    "Target Genes (comma-separated) ?", value="edd",
    help="Enter gene names (comma-separated) for spacer design."
)
pam_sequence = st.text_input(
    "PAM Sequence ?", value="CC",
    help="Specify PAM sequence (e.g., NGG for SpCas9)."
)
spacer_length = st.slider(
    "Spacer Length ?", 16, 50, 32,
    help="Slider to select spacer length."
)
direction = st.selectbox(
    "Spacer Direction ?",
    ["PAM-Upstream (PAM-Spacer)", "PAM-Downstream (Spacer-PAM)"],
    help="Choose whether PAM is upstream or downstream of spacer."
)
num_spacers = st.slider(
    "Number of Spacers to Predict per Gene ?", 1, 10, 3,
    help="How many spacers to predict per gene."
)
position_range = st.slider(
    "Search Position in Gene (% of Gene Length) ?", 0, 100, (10, 20),
    help="Percentage of gene length to search for spacers."
)

# Upstream design controls
design_upstream = st.checkbox(
    "Design Spacers Upstream of ORF ?",
    help="Enable to design spacers upstream of ORF start."
)
if design_upstream:
    upstream_window = st.slider(
        "Upstream Window for PAM (bases from ORF start) ?", 10, 500, (10, 100),
        help="Range upstream of ORF to search for PAM."
    )
else:
    upstream_window = None

# Flank sequences
left_flank = st.text_input(
    "Left Flank Sequence ?", value="aggtcTcaaaac",
    help="Sequence to prepend to each spacer."
)
right_flank = st.text_input(
    "Right Flank Sequence ?", value="gtttttGAGACCa",
    help="Sequence to append to each spacer."
)

# Function to extract ORF and upstream sequences
def get_orf_upstream_sequences(gb_content, gene_names, upstream_window, design_upstream):
    gene_seqs = {}
    upstream_seqs = {}
    handle = StringIO(gb_content)
    for record in SeqIO.parse(handle, "genbank"):
        for feat in record.features:
            if feat.type == "CDS" and "gene" in feat.qualifiers:
                name = feat.qualifiers["gene"][0]
                if name in gene_names:
                    # Extract ORF sequence
                    orf_seq = feat.location.extract(record).seq
                    gene_seqs[name] = orf_seq
                    # Extract upstream if requested
                    if design_upstream and upstream_window:
                        start = feat.location.start if feat.location.strand == 1 else feat.location.end
                        if feat.location.strand == 1:
                            u_start = max(0, start - upstream_window[1])
                            u_end = max(0, start - upstream_window[0])
                            up_seq = record.seq[u_start:u_end]
                        else:
                            u_start = min(len(record.seq), start + upstream_window[0])
                            u_end = min(len(record.seq), start + upstream_window[1])
                            up_seq = record.seq[u_start:u_end].reverse_complement()
                        upstream_seqs[name] = up_seq
    return gene_seqs, upstream_seqs

# Function to generate spacers with PAM
def generate_spacers_with_pam(gene_seqs, upstream_seqs, pam, length, direction, count, prange, design_up):
    out = {}
    source = upstream_seqs if design_up else gene_seqs
    for name, seq in source.items():
        seq = seq.upper()
        L = len(seq)
        start = int(prange[0]/100 * L)
        end = int(prange[1]/100 * L)
        found = []
        for i in range(start, min(end, L - length - len(pam)) + 1):
            if seq[i:i+len(pam)] == pam:
                if direction.startswith("PAM-Up"):  # PAM-Spacer
                    spac = seq[i+len(pam):i+len(pam)+length]
                else:  # Spacer-PAM
                    spac = seq[i-length:i] if i >= length else ""
                if spac:
                    found.append(spac)
                    if len(found) == count:
                        break
        out[name] = found
    return out

# Main action
if st.button("Generate Spacers ?", help="Generate spacers with current settings."):
    if not genbank_file:
        st.error("Please upload a GenBank file first.")
    else:
        gb_content = genbank_file.getvalue().decode()
        genes = [g.strip() for g in target_genes.split(',')]
        gseq, useq = get_orf_upstream_sequences(gb_content, genes, upstream_window, design_upstream)
        spc = generate_spacers_with_pam(gseq, useq, pam_sequence, spacer_length, direction, num_spacers, position_range, design_upstream)
        # Build export data with complement rows
        rows = []
        for gene, lst in spc.items():
            for idx, s in enumerate(lst, start=1):
                full = f"{left_flank}{s}{right_flank}"
                comp = str(Seq(full).reverse_complement())
                rows.append([f"{gene}_sp.{idx}", full])
                rows.append([f"{gene}_sp.{idx}_comp", comp])
        df = pd.DataFrame(rows, columns=["Name", "Sequence"])
        st.dataframe(df)
        # Export to Excel
        buf = BytesIO()
        with pd.ExcelWriter(buf, engine='openpyxl') as writer:
            df.to_excel(writer, index=False)
        buf.seek(0)
        st.download_button(
            "Download Results as Excel ?", buf, "spacers.xlsx",
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            help="Download spacers and complements as Excel"
        )
