# Repeated losses of self-fertility in _Saccharomyces cerevisiae_ evolution

## Description

This repository contains the data and scripts for the analyses presented in Vittorelli Nina, Gómez-Muñoz Cintia, Andriushchenko Irina, Ollivier Louis, Agier Nicolas, Delmas Stéphane, Corbeau Yann, Achaz Guillaume, Cosentino Lagomarsino Marco, Liti Gianni, Llorente Bertrand, Fischer Gilles, "Repeated losses of self-fertility shaped heterozygosity and polyploidy in yeast evolution".

## Installation instructions

First, download the repository with the download button, or clone it using the following command in your terminal.

```bash
git clone https://github.com/nina-vittorelli/Loss_SelfFertility_Scerevisiae.git
```

Second, install [conda](https://anaconda.org/anaconda/conda) if not already installed. Then open a terminal, and run the following commands in order to install the conda environment.

```bash
# go to the directory you downloaded or cloned
cd Loss_SelfFertility_Scerevisiae/

# create conda environment
conda env create -f environment/conda_env.yml -n loss_selfFertility_Scerevisiae_env

# activate conda environment
conda activate loss_selfFertility_Scerevisiae_env
```

Third, install [Rstudio](https://posit.co/downloads/) if not already installed. Then, install the `renv` R package (version 1.1.4) from a R terminal if not already installed. Finally, open the R project file `Loss_SelfFertility_Scerevisiae.Rproj`, and run the R notebooks inside this project. 

## Running the project

All scripts are in the folder `01_scripts`. They should be run in the order indicated by the numbers in the files' names. Bash and python scripts should be lauched from the terminal with the command indicated at the beginning of each file. R notebooks should be run within the R project `Loss_SelfFertility_Scerevisiae.Rproj`. 

Input data are in the `02_data` folder. A `README.md` file in this folder describes the origin of each file.

## License 

GNU GENERAL PUBLIC LICENSE v.3

## Author

Written by Nina Vittorelli

Contact: nina.vittorelli [at] sorbonne-universite.fr

