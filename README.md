# Diffusion coefficient predictor

Diffusion coefficients for several protein species are often required when modelling protein-protein interactions between several spatial compartments, however they are often difficult to find in the literature. As an alternative, diffusion coefficients can be estimated by fitting a mathematically formalised biochemical model to available data; however, an initial guess of the value of the diffusion coefficient is still required. Additionally, less time is required to fit models to data when the initial values are close to optimal. This program offers a quick and easy estimation of a likely range of diffusion coefficients based on the amino acid sequence of a protein.

## Getting started



### Requirements

The program requires Python 3 to successfully run. Additionally, several libraries are used by the program:

```
Bio
pandas
math
```
All required libraries are readily available for Pip installation.



### File types

The program accepts local .txt files with fasta sequences, or it can automatically download and process fasta files from the NCBI protein data base from a list of accession numbers.



## Using the pipeline

