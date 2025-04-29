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

# File uploads
genbank_file = st.file_uploader(
    "Upload a GenBank file (.gb or .gbff)", type=["gb", "gbff"],
    help="Upload the GenBank file containing your target gene's sequence."
)
gene_file = st.file_uploader(
    "Upload Gene List (Excel)", type=["xlsx", "xls"],
    help="Upload an Excel file with gene names in the first column. Overrides manual entry."
)

# Manual inputs
target_genes = st.text_input(
    "Target Genes (comma-separated)", value="edd",
    help="Enter gene names separated by commas. Ignored if Excel uploaded."
)
pam_sequence = st.text_input(
    "PAM Sequence", value="CC",
    help="Specify PAM sequence (e.g., NGG for SpCas9)."
)
spacer_length = st.slider(
    "Spacer Length", 16, 50, 32,
    help="Length of the spacer sequence."
)
direction = st.selectbox(
    "Spacer Direction",
    ["PAM-Upstream (PAM-Spacer)", "PAM-Downstream (Spacer-PAM)"],
    help="Orientation of spacer relative to PAM."
)
num_spacers = st.slider(
    "Number of Spacers per Gene", 1, 10, 3,
    help="How many spacers to predict per gene."
)
position_range = st.slider(
    "Search Position in Gene (% of Gene)", 0, 100, (10, 20),
    help="Percentage range within gene to search for spacers."
)

# Upstream design
design_upstream = st.checkbox(
    "Design Spacers Upstream of ORF",
    help="Design spacers in a window upstream of ORF start."
)
if design_upstream:
    upstream_window = st.slider(
        "Upstream Window (bases from ORF start)", 10, 500, (10, 100),
        help="Window upstream of ORF for PAM search."
    )
else:
    upstream_window = None

# Flank sequences
left_flank = st.text_input(
    "Left Flank Sequence", value="aggtcTcaaaac",
    help="Sequence to prepend to each spacer."
)
right_flank = st.text_input(
    "Right Flank Sequence", value="gtttttGAGACCa",
    help="Sequence to append to each spacer."
)

# Extract ORF/upstream sequences
def get_orf_upstream_sequences(gb_content, gene_list, upstream_window, design_up):
    gene_seqs, upstream_seqs = {}, {}
    handle = StringIO(gb_content)
    for rec in SeqIO.parse(handle, "genbank"):
        for feat in rec.features:
            if feat.type == "CDS" and "gene" in feat.qualifiers:
                name = feat.qualifiers["gene"][0]
                if name in gene_list:
                    seq = feat.location.extract(rec).seq
                    gene_seqs[name] = seq
                    if design_up and upstream_window:
                        start = feat.location.start if feat.location.strand == 1 else feat.location.end
                        if feat.location.strand == 1:
                            u0 = max(0, start - upstream_window[1])
                            u1 = max(0, start - upstream_window[0])
                            upstream_seqs[name] = rec.seq[u0:u1]
                        else:
                            u0 = min(len(rec.seq), start + upstream_window[0])
                            u1 = min(len(rec.seq), start + upstream_window[1])
                            upstream_seqs[name] = rec.seq[u0:u1].reverse_complement()
    return gene_seqs, upstream_seqs

# Generate spacers
def generate_spacers_with_pam(gene_seqs, upstream_seqs, pam, length, direction, count, prange, design_up):
    results = {}
    source = upstream_seqs if design_up else gene_seqs
    for name, seq in source.items():
        s = seq.upper()
        L = len(s)
        start, end = int(prange[0]/100*L), int(prange[1]/100*L)
        found = []
        for i in range(start, min(end, L-length-len(pam))+1):
            if s[i:i+len(pam)] == pam:
                spac = (s[i+len(pam):i+len(pam)+length] if direction.startswith("PAM-Up") else s[i-length:i] if i>=length else "")
                if spac:
                    found.append(spac)
                    if len(found) == count:
                        break
        results[name] = found
    return results

# Main action
if st.button("Generate Spacers", help="Generate spacers with current settings."):
    if not genbank_file:
        st.error("Please upload a GenBank file first.")
    else:
        # Read genes from Excel or manual
        if gene_file:
            try:
                dfg = pd.read_excel(gene_file)
                genes = dfg.iloc[:,0].dropna().astype(str).tolist()
            except Exception as e:
                st.error(f"Reading gene list failed: {e}")
                genes = []
        else:
            genes = [g.strip() for g in target_genes.split(',') if g.strip()]
        
        gb_text = genbank_file.getvalue().decode()
        gseq, useq = get_orf_upstream_sequences(gb_text, genes, upstream_window, design_upstream)
        spc = generate_spacers_with_pam(gseq, useq, pam_sequence, spacer_length, direction, num_spacers, position_range, design_upstream)
        
        # Build rows including complements
        rows = []
        for gene, lst in spc.items():
            for idx, seq in enumerate(lst, 1):
                full = f"{left_flank}{seq}{right_flank}"
                comp = str(Seq(full).reverse_complement())
                rows.append([f"{gene}_sp.{idx}", full])
                rows.append([f"{gene}_sp.{idx}_comp", comp])
        df = pd.DataFrame(rows, columns=["Name","Sequence"])
        st.dataframe(df)
        
        # Export
        buf = BytesIO()
        try:
            with pd.ExcelWriter(buf, engine='openpyxl') as writer:
                df.to_excel(writer, index=False)
            buf.seek(0)
            st.download_button("Download Results as Excel", buf, "spacers.xlsx", 
                                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        except Exception:
            csv_buf = BytesIO()
            df.to_csv(csv_buf, index=False)
            csv_buf.seek(0)
            st.download_button("Download Results as CSV", csv_buf, "spacers.csv", "text/csv")
