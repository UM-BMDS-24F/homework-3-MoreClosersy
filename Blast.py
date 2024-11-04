from Bio import SeqIO
from Bio.Blast import NCBIXML
import subprocess
import os
import io  

# File paths
human_file = "human.fa"
mouse_file = "mouse.fa"
mouse_db = "mouse_db"  

# Parsing human sequences
human_sequences = list(SeqIO.parse(human_file, "fasta"))

# Creating a BLAST database if it doesn't already exist
if not os.path.exists(mouse_db + ".psq"):
    subprocess.run(["makeblastdb", "-in", mouse_file, "-dbtype", "prot", "-out", mouse_db])

# Running BLAST for each human sequence
output_file = "blast_results.txt"
with open(output_file, "w") as f:
    for human_seq in human_sequences:
        with open("temp_human_query.fa", "w") as temp_query_file:
            SeqIO.write(human_seq, temp_query_file, "fasta")

        blast_command = [
            "blastp", 
            "-query", "temp_human_query.fa", 
            "-db", mouse_db, 
            "-evalue", "1e-5", 
            "-outfmt", "5", 
            "-num_alignments", "1",
            "-matrix", "BLOSUM62"
        ]
        
        # Executing BLAST
        result = subprocess.run(blast_command, capture_output=True, text=True)
        
        # Checking for errors
        if result.returncode != 0:
            print(f"Error running BLAST: {result.stderr}")
            continue

        # Parsing the BLAST output
        blast_record = NCBIXML.read(io.StringIO(result.stdout))

        # Processing results if alignments are found
        if blast_record.alignments:
            top_alignment = blast_record.alignments[0]
            top_hsp = top_alignment.hsps[0]
            
            # Writing the results to the output file
            f.write(f"Human Seq ID: {human_seq.id}\n")
            f.write(f"Mouse Seq ID: {top_alignment.hit_id}\n")  # Keeping the full ID
            f.write(f"Alignment:\n{top_hsp.sbjct}\n")
            f.write(f"E-value: {top_hsp.expect}\n")
            f.write(f"Bitscore: {top_hsp.bits}\n\n")
        
        # Remove the temporary file after each run
        os.remove("temp_human_query.fa")
            
# 1. I used BLASTp for the comparison since it’s designed for aligning protein sequences,
# making it suitable for analyzing similarities between mouse and human proteins.

# 2. BLOSUM62 was chosen as it is the default matrix for BLASTp and is particularly effective
# for proteins with moderate divergence, which is typical in comparisons between two mammalian species.

# 3. An E-value of 0.00001 ensures a high threshold, filtering out random matches while still capturing significant alignments.
# By setting the number of alignments to 1, we focus on the most homologous sequence, minimizing processing needs.
# Additionally, using the XML format (outfmt 5) simplifies parsing for IDs and alignment details.