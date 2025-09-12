#!/bin/zsh

##########################################################################
########## BUILD MULTIPLE SEQUENCE ALIGNMENT FOR HO CDS IN SCRAP #########
##########################################################################

# before running the script, activate the conda environment in your terminal with this command:
# conda activate loss_selfFertility_Scerevisiae_env
# To run this script:
# nohup zsh 01_scripts/08_scrap_HO_msa.sh 1>01_scripts/logs/log08a.txt 2>01_scripts/logs/log08b.txt &

# Paths to access the files
# input
PATH_SCRAP_ASSEMBLIES="02_data/01_scrap/nuclear_assemblies"
PHS2_SEQ="02_data/05_else/03_pHS2_sequence/HOcds_in_pHS2.fasta"
PATH_BED_FILES="03_output/05_HO_sequences_scrap/bedfiles"
# output
PATH_FASTA="03_output/05_HO_sequences_scrap/fastafiles"
CONCAT_FILE="$PATH_FASTA"/HOcds_scrap_nt.untrimmed.fasta
MERGED_FILE="$PATH_FASTA"/HOcds_scrap_pHS2_nt.untrimmed.fasta
ALIGNMENT="$PATH_FASTA"/HOcds_pHS2_scrap_nt.untrimmed.aligned.fasta
TRIMMED_ALIGNMENT="$PATH_FASTA"/HOcds_pHS2_scrap_nt.aligned.fasta
TRIMMED_UNALIGNED="$PATH_FASTA"/HOcds_pHS2_scrap_nt.fasta

# create the output directory 
mkdir -p "$PATH_FASTA"


##### I. EXTRACT SEQUENCES FROM BED #####

# loop over the bed files
for FILE in "$PATH_BED_FILES"/*_to_extract.bed; do

    # the path to the corresponding assembly
    ASSEMBLY="$PATH_SCRAP_ASSEMBLIES"/$(basename "${FILE%_to_extract.bed}")

    # the path to the folder for the output fasta sequences
        OUTPUT_FASTA="$PATH_FASTA"/$(basename "${FILE%_to_extract.bed}")_extracted.fasta

    # the bedtools command creates fasta from bed and assembly
    bedtools getfasta -fi "$ASSEMBLY" -bed "$FILE" -fo "$OUTPUT_FASTA" -s
        
done


##### II. CREATE A SINGLE FASTA FILE PER REGION #####

# remove the concatenated file present
rm "$CONCAT_FILE"

# create the concatenated file
touch "$CONCAT_FILE"

# add the file name at the begining of each sequence
for FILE in "$PATH_FASTA"/*.nuclear_genome.tidy.fa_extracted.fasta; do

    # strain + assembly informations
    FILENAME=$(basename "$FILE" | awk -F '.' '{print $1"."$3}')

    # add it to the sequences names
    awk -v filename="$FILENAME" '/^>/ {sub(/\([+-]\)$/, "", $1); gsub(":", "-", $1); sub(/^>/, ">"filename"-")} 1' "$FILE" >> "$CONCAT_FILE"
done 

# remove the individual files
rm "$PATH_FASTA"/*_extracted.fasta


##### III. ALIGN #####

# add the pHS2 sequences
cat "$PHS2_SEQ" "$CONCAT_FILE" > "$MERGED_FILE"

# multiple sequence alignment
mafft --preservecase --inputorder "$MERGED_FILE" | \
    # from multiple lines to single line for a sequence
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | \
    # remove the first line (empty)
    sed '1d' > "$ALIGNMENT"


##### IV. TRIM SEQUENCES #####

# Retrieve first sequence (HO_WT) 
PHS2_ALIGNED=$(awk '!/^>/ {print; exit}' "$ALIGNMENT")

# count
LEN_ALIGNMENT=$(echo "$PHS2_ALIGNED" | wc -c)
START_SITE=$(echo "$PHS2_ALIGNED" | sed 's/[^-].*//g' | wc -c)
END_SITE=$(echo "$PHS2_ALIGNED" | sed 's/-*$//' | wc -c)

echo "Length of the alignment before trimming: $LEN_ALIGNMENT"
echo "Trimmed alignment starts at site $START_SITE and ends and site $END_SITE."

# trim
awk -v Ni="$START_SITE" -v Nf="$END_SITE" '/^[^>]/ {print substr($0, Ni, Nf-Ni)} /^[>]/ {print}' "$ALIGNMENT" > "$TRIMMED_ALIGNMENT"

LEN_ALIGNMENT_TRIMMED=$(wc -L "$TRIMMED_ALIGNMENT")
echo "Length of the alignment after trimming: $LEN_ALIGNMENT_TRIMMED"

# Remove the gaps (unalign the sequences)
awk '/^>/ {print} !/^>/ {gsub("-", ""); print}' "$TRIMMED_ALIGNMENT" > "$TRIMMED_UNALIGNED"