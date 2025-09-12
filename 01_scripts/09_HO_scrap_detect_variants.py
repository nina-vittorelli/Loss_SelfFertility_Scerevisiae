##########################################################################
############# DETECT THE VARIANTS IN THE HO MSA OF THE SCRAP #############
##########################################################################

# before running the script, activate the conda environment in your terminal with this command:
# conda activate loss_selfFertility_Scerevisiae_env
# then run the script with the command:
# python 01_scripts/09_HO_scrap_detect_variants.py

### PACKAGES 
import re
import pandas as pd
import math
from Bio import AlignIO
from Bio import Seq



### PATHS TO FILES 
alignment_file = "03_output/05_HO_sequences_scrap/fastafiles/HOcds_pHS2_scrap_nt.aligned.fasta"

# for ambiguous bases
iupac_code = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"}}

### FUNCTIONS
def detect_variant(refseq, testseq):
    "Detect all variants in a sequence testseq aligned to a reference sequence refseq"

    variations = []

    # loop over positions
    i, s = 0, len(refseq)
    while i < s: 

        # no variant
        if refseq[i] == testseq[i]:
            i = i+1

        # deletion
        elif testseq[i] == "-":  

            # deletion size
            lendel = len(re.match(r"^-+", testseq[i:s]).group(0)) # count the number of "-"
            lendel_corrected = lendel - len(re.findall(r"-", refseq[i:i+lendel])) # remove the "-" in the ref

            # appends the variation to the list of variations
            variations.append([i + 1, "DEL", lendel_corrected, refseq[i:i+lendel].replace("-", ""),"*", 1])

            i = i + lendel

        # insertion
        elif refseq[i] == "-":  

            # insertion size
            lenins = len(re.match(r"^-+", refseq[i:s]).group(0)) # count the number of "-" in the ref
            lenins_corrected = lenins - len(re.findall(r"-", testseq[i:i+lenins])) # remove the "-" of the tested sequence

            # insertion sequence
            seqins = testseq[i:i+lenins].replace("-", "")

            # appends the variation to the list of variations
            variations.append([i + 1, "INS", lenins_corrected, "*", seqins, 1])

            i = i + lenins

        # SNP
        else: 

            if testseq[i] in iupac_code:
                possible_bases = iupac_code[testseq[i]]

                if refseq[i] in possible_bases:
                    # Partial mismatch — find the base that does NOT match the reference
                    alt_base = (possible_bases - {refseq[i]}).pop()
                    variations.append([i + 1, "SNP", 1, refseq[i], alt_base, 0.5])

                else: 
                    for base in possible_bases: 
                        variations.append([i + 1, "SNP", 1, refseq[i], base, 0.5])
            else: 
                variations.append([i + 1, "SNP", 1, refseq[i], testseq[i], 1])

            i = i + 1

    return variations

def build_variants_table(alignment, ref, quiet): 
    "Build a dataframe with all the variants from the ref sequence in an alignment."

    # create an empty dataframe
    df = pd.DataFrame(columns=["pos", "type", "size", "ref_nt", "alt_nt"])
    
    refseq = str(alignment[id == ref].seq).upper()

    # loop over the sequences
    for sequence in alignment: 
        
        # skip reference
        if sequence.id == ref: 
            continue
        
        else: 
            # detect the variant of the sequence 
            testid, testseq = sequence.id, str(sequence.seq).upper()
            variants = detect_variant(refseq, testseq) 

            # create a dataframe
            variants_df = pd.DataFrame(variants, columns=["pos", "type", "size", "ref_nt", "alt_nt", testid])

            # print info
            nVariants = variants_df.shape[0]
            nSNPs = variants_df.loc[variants_df.type == "SNP"].shape[0]
            nINSs = variants_df.loc[variants_df.type == "INS"].shape[0]
            nDELs = variants_df.loc[variants_df.type == "DEL"].shape[0]
            
        
            if quiet == False: 
                print(testid + ":\t" + str(nVariants) + " variants, including \t" + str(nSNPs) + " SNPs, \t" + str(nINSs) + " insertions and \t" + str(nDELs) + " deletions" )
            
            # merge the results to the previous dataframe
            df = pd.merge(df, variants_df, on=["pos", "type", "size","ref_nt", "alt_nt"], how="outer")

    df = df.fillna(0)

    # indicate the position of the refseq (in case of insertions)
    df["pos_ref"] = [i - len(re.findall(r"-", refseq[0:i])) for i in df["pos"]]

    # print info
    nVariants = df.shape[0]
    nSequences = df.shape[1] - 4
    nSNPs = df.loc[df.type == "SNP"].shape[0]
    nINSs = df.loc[df.type == "INS"].shape[0]
    nDELs = df.loc[df.type == "DEL"].shape[0]

    print("Total number of sequences: " + str(nSequences))
    print("Total:\t" + str(nVariants) + " variants, including \t" + str(nSNPs) + " SNPs, \t" + str(nINSs) + " insertions and \t" + str(nDELs) + " deletions" )

    return df


def variant_effect_cds(alignment, ref="default", quiet=False): 

    # set up the reference sequence
    if ref == "default": 
        ref = alignment[0].id 
    print("Reference: " + str(ref)) 

    refseq = str(alignment[id == ref].seq).upper()
    refseq_no_gap = re.sub(r"-", "", refseq)

    if re.sub(r"-", "", refseq).startswith("ATG") == False: 
        print("ERROR: refseq should start with a start codon ATG")

    if (re.sub(r"-", "", refseq).endswith("TAA") or re.sub(r"-", "", refseq).endswith("TAG")  or re.sub(r"-", "", refseq).endswith("TGA")) == False: 
        print("ERROR: refseq should start with a start codon ATG")

    print("length alignment: " + str(len(refseq)))
    print("nb of - in ref: " + str(len(re.findall(r"-", refseq))))
    print("len of refseq (nt): " + str(len(refseq) - len(re.findall(r"-", refseq))))
    
    if (len(refseq) - len(re.findall(r"-", refseq))) % 3 == 0:
        refProtLen = int((len(refseq) - len(re.findall(r"-", refseq))) / 3)
        print("Length of the ref protein (aa): " + str(refProtLen))
    
    else: 
        print("ERROR: refseq is not a coding sequence (length is not divisible by 3)")
        
    # search for variants
    variants_df = build_variants_table(alignment, ref, quiet)

    # giving position of the variant relative to the coding feature
    variants_df["pos_codon"] = [math.ceil(i/3) for i in variants_df["pos_ref"]]
    variants_df["pos_in_codon"] = (variants_df["pos_ref"] % 3).replace(0,3)
    variants_df["ref_codon"] = [refseq_no_gap[(i-1) * 3 : i * 3] for i in variants_df["pos_codon"]]
    variants_df["ref_amino_acid"] = [Seq.Seq(codon).translate() for codon in variants_df["ref_codon"]]

    # loss start codon if there is a mutation in the 3 first nucleotides
    variants_df["loss_start_codon"] = (variants_df["pos_codon"] == 1)

    # in frame indels
    variants_df["in_frame_indel_between_codons"] = (variants_df["type"].isin(["DEL", "INS"])) & (variants_df["size"] % 3 == 0) & (variants_df["pos_in_codon"] == 1) 

    # in frame indels but within a codon
    variants_df["in_frame_indel_within_codon"] = (variants_df["type"].isin(["DEL", "INS"])) & (variants_df["size"] % 3 == 0) & (variants_df["pos_in_codon"] != 1) 

    # frameshift indels
    variants_df["frameshift"] = (variants_df["type"].isin(["DEL", "INS"])) & (variants_df["size"] % 3 != 0)

    # SNPs
    variants_df["alt_codon"] = variants_df.apply(
        lambda row: row["ref_codon"][:row["pos_in_codon"] - 1] + row["alt_nt"] + row["ref_codon"][row["pos_in_codon"]:] if row["type"] == "SNP" else "", 
        axis=1)
    variants_df["alt_amino_acid"] = [Seq.Seq(codon).translate() for codon in variants_df["alt_codon"]]

    variants_df["synonymous_SNP"] = (variants_df["type"] == "SNP") & (variants_df["loss_start_codon"] == False) & (variants_df["alt_amino_acid"] == variants_df["ref_amino_acid"])
    variants_df["nonsense_SNP"] = (variants_df["type"] == "SNP") & (variants_df["loss_start_codon"] == False) & (variants_df["alt_amino_acid"] != variants_df["ref_amino_acid"]) & (variants_df["alt_amino_acid"] == "*")
    variants_df["missense_SNP"] = (variants_df["type"] == "SNP") & (variants_df["loss_start_codon"] == False) & (variants_df["alt_amino_acid"] != variants_df["ref_amino_acid"]) & (variants_df["alt_amino_acid"] != "*")

    return variants_df


### MAIN 
# import the msa
alignment = AlignIO.read(alignment_file, "fasta")

# finds the variants in the alignement
table_variants_and_effect = variant_effect_cds(alignment, ref="HOgene_in_pHS2", quiet=False)

# write the table to csv
table_variants_and_effect.to_csv("03_output/05_HO_sequences_scrap/HO_variants_pHS2_vs_scrap.csv")