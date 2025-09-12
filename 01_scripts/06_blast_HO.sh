#!/bin/bash

##########################################################################
############# BLAST THE PHS2 SEQUENCE IN THE SCRAP ASSEMBLIES ############
##########################################################################

# before running the script, activate the conda environment in your terminal with this command:
# conda activate loss_selfFertility_Scerevisiae_env
# To run this script:
# nohup bash 01_scripts/06_blast_HO.sh  1>01_scripts/logs/log06a.txt 2>01_scripts/logs/log06b.txt &

# Paths to access the files

# Inputs
PATH_SCRAP_ASSEMBLIES="02_data/01_scrap/nuclear_assemblies"
HO_PATH="02_data/05_else/03_pHS2_sequence/HOcds_in_pHS2.fasta"

# Outputs
OUTPUT_FOLDER="03_output/05_HO_sequences_scrap/"
MERGED_OUTPUT="$OUTPUT_FOLDER"blast_merged.csv

# create the output directory 
mkdir -p "$OUTPUT_FOLDER"

# Uncompress the files if needed
if ls "$PATH_SCRAP_ASSEMBLIES"/*.fa.gz >/dev/null 2>&1; then
	for ASSEMBLY_FA_GZ in "$PATH_SCRAP_ASSEMBLIES"/*.fa.gz; do
		echo "Uncompress $(basename $ASSEMBLY_FA_GZ)"
		gunzip $ASSEMBLY_FA_GZ
	done
fi 

# Blast the sequence in pHS2 in the ScRAP assemblies
echo -e "assembly_file qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" > "$MERGED_OUTPUT"

for ASSEMBLY_FA in "$PATH_SCRAP_ASSEMBLIES"/*.fa; do

	# Create a database for blast if needed
	if [ ! -f "$ASSEMBLY_FA".ndb ]; then
		echo "Create the database for $(basename "$ASSEMBLY_FA")"
		makeblastdb -in "$ASSEMBLY_FA" -dbtype nucl -out "$ASSEMBLY_FA"
	fi
	
	echo "Blast in $(basename "$ASSEMBLY_FA")"
	# Blast
	blastn -query "$HO_PATH" \
		-db "$ASSEMBLY_FA" \
		-out "$OUTPUT_FOLDER"blast_$(basename "$ASSEMBLY_FA").csv \
		-num_threads 4 \
		-outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore"

	# create a single file with all the results
	while IFS=$'\t' read -r -a line; do
        echo -e "$(basename "$ASSEMBLY_FA") ${line[*]}" >> "$MERGED_OUTPUT"
    done < "$OUTPUT_FOLDER"blast_$(basename "$ASSEMBLY_FA").csv

	# remove the per-assembly files
	rm "$OUTPUT_FOLDER"blast_$(basename "$ASSEMBLY_FA").csv
done
