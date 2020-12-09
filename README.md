# Diffusion coefficient predictor

Diffusion coefficients for several protein species are often required when modelling protein-protein interactions between several spatial compartments, however they are often difficult to find in the literature. As an alternative, diffusion coefficients can be estimated by fitting a mathematically formalised biochemical model to available data; however, an initial guess of the value of the diffusion coefficient is still required. Additionally, less time is required to fit models to data when the initial values are close to optimal. This program offers a quick and easy estimation of a likely range of diffusion coefficients based on the amino acid sequence of a protein.

## Getting started



### Requirements

The program requires Python 3 to successfully run. Additionally, several libraries are used by the program:

```
Bio
pandas
math
sys
```
All required libraries are readily available for Pip installation.



### File types

The program accepts local .txt files with fasta sequences, or it can automatically download and process fasta files from the NCBI protein data base from a list of accession numbers.



## Using the program



### Testing and background

After installation the program can be tested by running the TESTME.py file from a command line environment as follows:

```
$ python TESTME.py 'my_email@address.domain'
```

If the program successfully runs, it should output the following text:

```
The estimated range of diffusion coefficients for bovine serum albumin is 0.66-16.54 um^2/s
The experimentally determined diffusion coefficient for BSA is 5.98 um^2/s

The estimated range of diffusion coefficients for chicken lysozyme is 1.03-25.81 um^2/s
The experimentally determined diffusion coefficient for chicken lysozyme is 11.4 um^2/s

The estimated range of diffusion coefficients for yeast phosphoglycerate kinase is 0.76-18.9 um^2/s
The experimentally determined diffusion coefficient for YPGK is 6.38 um^2/s
```

The script relies on biopython to access and use NCBI's Entrez tools, which require users to submit an email address in case of overloading servers. Biopython has built-in safeguards against overloading servers, however failure to provide a valid email address may result in the banning of the offending IP address in case any problems arise. For more informations, see NCBI's [Entrez Programming Utilities resources](https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.chapter2_table1).

The test script probes NCBI's protein database for fasta sequences associated with three proteins with experimentally determined diffusion coefficients (from Tyn and Gusek, 1990). The program then calculates the molecular weight of each protein based on its amino acid sequence, which is then used to estimate its hydrodynamic radius. The function for estimating its hydrodynamic radius was developed in the jupyter notebook `Analyse-data_protein-weight-dimension.ipynb` with data from Tyn and Gusek (1990). Using the hydrodynamic radius, a range of theoretical diffusion coefficients at 37&deg;C is calculated using the [Stokes-Einstein equation](https://en.wikipedia.org/wiki/Einstein_relation_(kinetic_theory). The range is estimated by varying the viscosity of the cytoplasm of the cell between 2 and 50 cps (Obodovskiy, 2019), all other variables remain constant.



### Using the program with accession numbers


