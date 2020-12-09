# DCE: Diffusion coefficient estimation

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
>The estimated range of diffusion coefficients for bovine serum albumin is 0.66-16.54 um^2/s
>The experimentally determined diffusion coefficient for BSA is 5.98 um^2/s
>
>The estimated range of diffusion coefficients for chicken lysozyme is 1.03-25.81 um^2/s
>The experimentally determined diffusion coefficient for chicken lysozyme is 11.4 um^2/s
>
>The estimated range of diffusion coefficients for yeast phosphoglycerate kinase is 0.76-18.9 um^2/s
>The experimentally determined diffusion coefficient for YPGK is 6.38 um^2/s
```

The script relies on biopython to access and use NCBI's Entrez tools, which require users to submit an email address in case of overloading servers. Biopython has built-in safeguards against overloading servers, however failure to provide a valid email address may result in the banning of the offending IP address in case any problems arise. For more informations, see NCBI's [Entrez Programming Utilities resources](https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.chapter2_table1).

The test script probes NCBI's protein database for fasta sequences associated with three proteins with experimentally determined diffusion coefficients (from Tyn and Gusek, 1990). The program then calculates the molecular weight of each protein based on its amino acid sequence, which is then used to estimate its hydrodynamic radius. The function for estimating its hydrodynamic radius was developed in the jupyter notebook `Analyse-data_protein-weight-dimension.ipynb` with data from Tyn and Gusek (1990). Using the hydrodynamic radius, a range of theoretical diffusion coefficients at 37&deg;C is calculated using the [Stokes-Einstein equation](https://en.wikipedia.org/wiki/Einstein_relation_(kinetic_theory). The range is estimated by varying the viscosity of the cytoplasm of the cell between 2 and 50 cps (Obodovskiy, 2019), all other variables remain constant.



### Using the program with accession numbers

The program allows the importation and analysis of several proteins directly from a single function based on a list of accession numbers and outpts a pandas dataframe with the result. While two modules are needed, only one is required to be directly imported. See example below on how to use the program for estimation of diffusion coefficients of three proteins.

```
import diffusion_estimate as de

#Example list of accession numbers for BSA, chicken lysozyme and yeast phosphoglycerate kinase
example = = ['CAA76847.1', 'ACL81750.1',  'CAA42329.2']

data = de.DiffusionCoefficient('my_email@address.domain', example')

data

	ID	kDa	R, nm	D max, um^2/s	D min, um^2/s
0	CAA76847.1	80.2424	6.86669	16.5415	0.661662
1	ACL81750.1	18.7935	4.40146	25.8064	1.03225
2	CAA42329.2	52.2154	6.00911	18.9023	0.75609
```


### Using the program with a local text file of fasta sequences

The program can also be used with local .txt files of fasta sequences. Two example files with a single (test.fasta) and multiple (test2.fasta) fasta sequences are provided in the main branch. For local files, no email is required.

```
import diffusion_estimate as de

data = de.DiffusionCoefficient(my_email = None, acc_numbers = None, download = False, file = 'test2.fasta')

data

	ID	kDa	R, nm	D max, um^2/s	D min, um^2/s
0	WP_012282380.1	78.6474	6.82383	16.6455	0.665818
1	5V8K_A	77.6455	6.7966	16.7121	0.668485
2	WP_161256521.1	78.426	6.81783	16.6601	0.666404
3	WP_161261372.1	78.5532	6.82128	16.6517	0.666067
4	WP_155475590.1	78.8849	6.83025	16.6298	0.665193
5	AAA02821.1	78.5929	6.82235	16.6491	0.665962
6	WP_151619721.1	78.6343	6.82347	16.6463	0.665853
7	QGG48731.1	78.8487	6.82927	16.6322	0.665288
8	WP_162008027.1	75.3339	6.73288	16.8703	0.674812
9	WP_131918022.1	78.7661	6.82704	16.6376	0.665505
10	WP_155475454.1	59.5632	6.25896	18.1477	0.725908
```


### Limitations

Currently the program does not support fasta files with gaps in the sequence, and is therefore best suited to analysis of curated protein sequences.



## References

Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25, 1422-1423

McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*. 445.

Obodovskiy, I. (2019). Basics of Biology. *Radiation*: 429-445.

The Pandas Development Team (2020). pandas-dev/pandas: Pandas. [10.5281/zenodo.3509134](10.5281/zenodo.350913)

Tyn, M. T. and T. W. Gusek (1990). "Prediction of diffusion coefficients of proteins." *Biotechnology and Bioengineering* 35(4): 327-338.

Van Rossum, G., 2020. The Python Library Reference, release 3.8.2, Python Software Foundation.
