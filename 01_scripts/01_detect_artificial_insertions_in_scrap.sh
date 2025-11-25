#!/bin/bash

##########################################################################
############ BLAST THE MARKER SEQUENCES IN THE SCRAP ASSEMBLIES ##########
##########################################################################

# before running the script, activate the conda environment in your terminal with this command:
# conda activate loss_selfFertility_Scerevisiae_env
# To run this script:
# nohup bash 01_scripts/01_detect_artificial_insertions_in_scrap.sh 1>logs/log01a.txt 2>logs/log01b.txt &

# Paths
PATH_SCRAP_ASSEMBLIES="02_data/01_scrap/nuclear_assemblies"
PATH_MARKERS="02_data/05_else/01_artificial_markers"
PATH_OUTPUT="03_output/01_cleaned_panels/01_scrap_assemblies_with_markers"

# create the output directory 
mkdir -p "$PATH_OUTPUT"

# uncompress the files of the ScRAP genome assemblies if needed
if ls "$PATH_SCRAP_ASSEMBLIES"/*.fa.gz >/dev/null 2>&1; then
	for ASSEMBLY_FA_GZ in "$PATH_SCRAP_ASSEMBLIES"/*.fa.gz; do
		echo "Uncompress $(basename $ASSEMBLY_FA_GZ)"
		gunzip $ASSEMBLY_FA_GZ
	done
fi 

# loop over the genetic markers in the folder PATH_MARKERS
for MARKER_FASTA in "$PATH_MARKERS"/*.fasta; do 
	MARKER="$(basename "$MARKER_FASTA")"
	MARKER="${MARKER%%_nt.fasta}"

	echo "Blasting $MARKER in the assemblies"

	# create a csv file that will store the results
	MERGED_OUTPUT="$PATH_OUTPUT"/_"$MARKER"_scrap_BLAST.csv
	echo -e "assembly_file qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" > "$MERGED_OUTPUT"

	# loop over the genome assemblies of blast
	for ASSEMBLY_FA in "$PATH_SCRAP_ASSEMBLIES"/*.fa; do

		# create a database for blast if needed
		if [ ! -f "$ASSEMBLY_FA".ndb ]; then
			echo "Create the database for $(basename "$ASSEMBLY_FA")"
			makeblastdb -in "$ASSEMBLY_FA" -dbtype nucl -out "$ASSEMBLY_FA"
		fi
		
		# blast the marker sequence - create a file per strain x marker
		blastn -query "$MARKER_FASTA" \
			-db "$ASSEMBLY_FA" -out "$PATH_OUTPUT"/"$MARKER"_$(basename "$ASSEMBLY_FA").csv \
			-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

		# add the information to the marker csv file. 
		while IFS=$'\t' read -r -a line; do
        	echo -e "$(basename "$ASSEMBLY_FA") ${line[*]}" >> "$MERGED_OUTPUT"
    	done < "$PATH_OUTPUT"/"$MARKER"_$(basename "$ASSEMBLY_FA").csv

        rm "$PATH_OUTPUT"/"$MARKER"_$(basename "$ASSEMBLY_FA").csv

    done

done
