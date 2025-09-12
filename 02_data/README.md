# Origins of the data

## 01. ScRAP

### 01. Nuclear assemblies (folder nuclear_assemblies)
Downloaded from https://www.evomicslab.org/db/ScRAPdb/download/ on the 21-10-2024.

### 02. Metadata (file ScRAP_db.xlsx)
Metadata downloaded from https://www.nature.com/articles/s41588-023-01459-y#Sec40 (Supplementary data of: O’Donnell S, Yue J-X, Saada OA, Agier N, Caradec C, Cokelaer T, et al. Telomere-to-telomere assemblies of 142 strains characterize the genome structural landscape in Saccharomyces cerevisiae. Nat Genet. 2023;55: 1390–1399. doi:10.1038/s41588-023-01459-y)

### 03 Phylogenetic tree (file scrap_tree.tre)
Phylogenetic tree of the ScRAP strains from the original paper (O’Donnell S, Yue J-X, Saada OA, Agier N, Caradec C, Cokelaer T, et al. Telomere-to-telomere assemblies of 142 strains characterize the genome structural landscape in Saccharomyces cerevisiae. Nat Genet. 2023;55: 1390–1399. doi:10.1038/s41588-023-01459-y)


## 02. 1011 Genomes Panel (file SupplTables_1011Genomes.xls) 
Metadata downloaded from https://www.nature.com/articles/s41586-018-0030-5#Sec40 (supplementary data of Peter J, De Chiara M, Friedrich A, Yue J-X, Pflieger D, Bergström A, et al. Genome evolution across 1,011 Saccharomyces cerevisiae isolates. Nature. 2018;556: 339–344. doi:10.1038/s41586-018-0030-5)


## 03. 3,034 Genomes Panel

### 01. Metadata (file Table_S1_G3-2024-405400.xlsx)
Metadata downloaded from https://academic.oup.com/g3journal/article/14/12/jkae245/7904545#supplementary-data (supplementary data from  Loegler V, Friedrich A, Schacherer J. Overview of the Saccharomyces cerevisiae population structure through the lens of 3,034 genomes. G3 GenesGenomesGenetics. 2024;14: jkae245. doi:10.1093/g3journal/jkae245)

### 02. Ecological Niches
For each strain, an ecological niche was attributed based on the keywords of the "Isolation" description (see keywords in the Materials and Methods section of the paper).

### 03. vcf filtered for HO (folder vcf_filter_HO)
- The VCF file of the HO locus (chromosome 4, positions 46271 to 48031) from the 3,034GP (Loegler et al. 2024) was provided by Jing Hou (Université de Strasbourg, CNRS, GMGM UMR 7156, Strasbourg, France). 
- the reference genome of S288C used for the VCF is in fna format in the same folder, from the Zenodo repository of Loegler et al. 2024 (https://doi.org/10.5281/zenodo.12580561). 


## 04. Experimental data (file experiments.xlsx)
Results from the experimental screen. Description of the sheets: 
- scrapOI: mating type, sporulation and ploidy of ScRAP non-monosporic isolates
- scrapMI: mating type and sporulation of ScRAP monosporic isolates
- scrapMI_OI: mating type and sporulation of the diploid parent of ScRAP monosporic isolates
- newMI_from_scrapOI: mating type and sporulation newly-created monoporic isolates of ScRAP non-monosporic isolates
- newMI_from_scrapMI_OI: mating type and sporulation newly-created monoporic isolates of from the parents of ScRAP monosporic isolates
- rescue_pHS2: result of the molecular complementation assay of heterothallic isolates with pHS2.


## 05. Else

### 01. Artificial markers
Sequences of common markers used to artificially knockout genes like HO. Sequences were downloaded from the NCBI, the identifiants and positions can be find in the id of each sequence. 

### 02. Refseqs for MAT, HML and HMR
Sequences extracted from the SGD website (https://www.yeastgenome.org/). Id and position in the SGDref genome was kept in the id of each sequence. 
X region and Z1 regions are common to MATa, MATalpha, HML and HMR.
Y region can be Ya (in MATa and HMRa) or Yalpha (MATalpha and HMLalpha). Y region is flanked by X and Z1 regions, which are common to MATa, MATalpha, HML and HMR.
VBA3 (YCL069W) and CHA1 (YCL064C) suround HMLalpha.
PHO87 (YCR037C) and TAF2 (YCR042C) suround MAT (a or alpha).
CDC50 (YCR094W) and GIT1 (YCR098C) suround HMRa.

### 03. pHS2 sequence
pHS2.fasta: Fasta file corresponding to the pHS2 plasmid (Addgene plasmid #81037 ; http://n2t.net/addgene:81037 ; RRID:Addgene_81037) that was used during the complementation test. 
HOcds_in_pHS2.fasta: coding sequence of HO in the pHS2 plasmid.


## 06 Figures
Figures built manually (schematics, annotated HO network).
