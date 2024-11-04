from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML
import os

# Paths to the human and mouse fasta files
human_fasta_path = "human.fa"
mouse_fasta_path = "mouse.fa"
output_file = "blast_results.txt"

# Prepare the BLAST database from mouse.fa
makeblastdb_cmd = f"makeblastdb -in {mouse_fasta_path} -dbtype nucl"
os.system(makeblastdb_cmd)

# Open the output file for writing results
with open(output_file, "w") as output:
    # Iterate through each sequence in the human file
    for record in SeqIO.parse(human_fasta_path, "fasta"):
        # Write the individual human sequence to a temporary file
        with open("temp_human.fa", "w") as temp_file:
            SeqIO.write(record, temp_file, "fasta")

        # Run the BLAST search
        blastn_cline = NcbiblastnCommandline(
            query="temp_human.fa",
            db=mouse_fasta_path,
            evalue=1e-5,
            outfmt=5,  # XML format for easy parsing
            out="temp_blast.xml"
        )
        stdout, stderr = blastn_cline()

        # Parse the BLAST XML output
        with open("temp_blast.xml") as blast_result:
            blast_records = NCBIXML.parse(blast_result)
            for blast_record in blast_records:
                if blast_record.alignments:
                    # Process the top alignment
                    top_hit = blast_record.alignments[0]
                    hsp = top_hit.hsps[0]
                    output.write(
                        f"Human ID: {record.id}\n"
                        f"Mouse ID: {top_hit.hit_id}\n"
                        f"Alignment:\n{hsp.sbjct}\n{hsp.query}\n\n"
                        f"E-value: {hsp.expect}\n"
                        f"Bitscore: {hsp.bits}\n\n"
                    )
                else:
                    output.write(f"Human ID: {record.id}\nNo homolog found.\n\n")

# Cleanup temporary files
os.remove("temp_human.fa")
os.remove("temp_blast.xml")
