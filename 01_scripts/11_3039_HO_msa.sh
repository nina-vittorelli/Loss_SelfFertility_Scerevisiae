#!/bin/zsh

##########################################################################
########## BUILD MULTIPLE SEQUENCE ALIGNMENT FOR HO IN THE 3039 ##########
##########################################################################

# before running the script, activate the conda environment in your terminal with this command:
# conda activate loss_selfFertility_Scerevisiae_env
# To run this script:
# nohup zsh 01_scripts/11_3039_HO_msa.sh 1>logs/log11a.txt 2>logs/log11b.txt &

# Input
PATH_REF_GENOME="02_data/03_3039GP/vcf_filtered_HO/Sace_S288c_reference_FullMatrixID.fna"
PATH_VCF="02_data/03_3039GP/vcf_filtered_HO/full3039Matrix.chr4_46271-48031_filtered.vcf.gz"
STRAINS_TO_REMOVE="03_output/01_cleaned_panels/strains_to_exclude_3034GP.txt"
REF_SEQ_PHS2="02_data/05_else/03_pHS2_sequence/HOcds_in_pHS2.fasta"

# Output
OUTPUT_FOLDER="03_output/06_HO_sequences_3039"
SAMPLES_ID="$OUTPUT_FOLDER"/samples.txt
FASTA_FROM_VCF="$OUTPUT_FOLDER"/HOcds_3039.fasta
REVERSE_COMPLEMENTED_FASTA="$OUTPUT_FOLDER"/HOcds_3039.rc.fasta
FASTA_FILTERED="$OUTPUT_FOLDER"/HOcds_filtered.fasta
FASTA_WITH_REF="$OUTPUT_FOLDER"/HOcds_filtered_with_ref.fasta
FINAL_FASTA="$OUTPUT_FOLDER"/HOcds_filtered_with_ref.aligned.fasta

# create the output directory 
mkdir -p "$OUTPUT_FOLDER"

# Position of HO
REGION="chromosome4:46271-48031"

# create the index and dictionnary for the reference assembly
samtools faidx "$PATH_REF_GENOME"
gatk CreateSequenceDictionary -R "$PATH_REF_GENOME"

# create the index for the vcf + store the name of the strains
gatk IndexFeatureFile -I "$PATH_VCF"
bcftools query -l "$PATH_VCF" > "$SAMPLES_ID"

touch "$FASTA_FROM_VCF"

while read SAMPLE; do

    echo "Processing sample: $SAMPLE"

    # reconstruct the sequence
    samtools faidx "$PATH_REF_GENOME" "$REGION" | \
    bcftools consensus "$PATH_VCF" -s "$SAMPLE" -p "$SAMPLE"_YDL227C_HO >> "$FASTA_FROM_VCF"

done < "$SAMPLES_ID"

seqtk seq -r "$FASTA_FROM_VCF" > "$REVERSE_COMPLEMENTED_FASTA"

SEQ_COUNT=$(grep -c "^>" "$REVERSE_COMPLEMENTED_FASTA")
echo "Number of sequences before filtering: $SEQ_COUNT"

# remove strains with artificial deletion of HO
PATTERN=$(sed 's/$/_/' "$STRAINS_TO_REMOVE" | paste -sd',' -)
seqkit grep -v -r -p "$PATTERN" "$REVERSE_COMPLEMENTED_FASTA" > "$FASTA_FILTERED"

SEQ_COUNT=$(grep -c "^>" "$FASTA_FILTERED")
echo "Number of sequences after filtering: $SEQ_COUNT"

cat  "$REF_SEQ_PHS2" "$FASTA_FILTERED" > "$FASTA_WITH_REF"

mafft --preservecase --inputorder "$FASTA_WITH_REF" | \
        # from multiple lines to single line for a sequence
        awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | \
        # remove the first line (empty)
        sed '1d' > "$FINAL_FASTA"