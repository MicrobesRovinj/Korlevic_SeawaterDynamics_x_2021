## Temporal variation in the prokaryotic community of a nearshore marine environment
This is the repository for the manuscript "Temporal variation in the prokaryotic community of a nearshore marine environment" written by Marino Korlević, Marsej Markovski, Gerhard J. Herndl and Mirjana Najdek. The raw sequencing data with the exception of negative controls have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession numbers SAMEA6648771 – SAMEA6648788, SAMEA6648824, SAMEA6648825, SAMEA8117500 – SAMEA8117516. Negative control samples are part of this repository and located in data/raw/. To be able to reproduce the results the mothur compatible SILVA reference file (Release 138, available under a [CC-BY license](https://www.arb-silva.de/silva-license-information/)) must be created according to the instruction given on [the mothur blog](https://mothur.org/blog/2020/SILVA-v138-reference-files/) and in the Makefile. This README file contains an overview of the repository structure, information on software dependencies and instructions how to reproduce and rerun the analysis.

### Overview

	project
	|- README                       # the top level description of content (this doc)
	|- LICENSE                      # the license for this project
	|
	|- submission/                  # files necessary for manuscript or supplementary information rendering, e.g. executable Rmarkdown
	| |- manuscript.Rmd             # executable Rmarkdown for the manuscript of this study
	| |- manuscript.md              # Markdown (GitHub) version of the manuscript.Rmd file
	| |- manuscript.tex             # TeX version of manuscript.Rmd file
	| |- manuscript.pdf             # PDF version of manuscript.Rmd file
	| |- manuscript.aux             # auxiliary file of the manuscript.tex file, used for cross-referencing
	| |- header.tex                 # LaTeX header file to format the PDF version of manuscript
	| |- supplementary.Rmd          # executable Rmarkdown for the supplementary information of this study
	| |- supplementary.md           # Markdown (GitHub) version of the supplementary.Rmd file
	| |- supplementary.tex          # TeX version of supplementary.Rmd file
	| |- supplementary.pdf          # PDF version of supplementary.Rmd file
	| |- supplementary.aux          # auxiliary file of the supplementary.tex file, used for cross-referencing
	| |- header_supplementary.tex   # LaTeX header file to format the PDF version of supplementary information
	| |- references.bib             # BibTeX formatted references
	| +- citation_style.csl         # csl file to format references
	|
	|- data                         # raw and primary data, are not changed once created
	| |- references/                # reference files to be used in analysis
	| |- raw/                       # raw data, will not be altered
	| |- mothur/                    # mothur processed data
	|
	|- code/                        # any programmatic code
	|
	|- results                      # all output from workflows and analyses
	| |- figures/                   # graphs designated for manuscript or supplementary information figures
	| +- numerical/                 # results of the statistics or other numerical results for manuscript or supplementary information
	|
	|-.gitignore                    # gitinore file for this study
	|-.Rprofile                     # Rprofile file containing information on which R libraries to load,
	|                               # rendering options for knitr and Rmarkdown, functions etc.
	+- Makefile                     # executable Makefile for this study

### How to regenerate this repository

#### Dependencies
* GNU Bash (v. 4.2.46(2))
* GNU Make (v. 3.82); should be located in the user's PATH
* mothur (v. 1.43.0)
* R (v. 3.6.0); should be located in the user's PATH
* R packages:
  * `stats (v. 3.6.0)`
  * `knitr (v. 1.25)`
  * `rmarkdown (v. 1.16)`
  * `tinytex (v. 0.16)`
  * `tidyverse (v. 1.2.1)`
  * `vegan (v. 2.5.6)`
  * `RColorBrewer (v. 1.1.2)`
  * `kableExtra (v. 1.1.0)`
  * `cowplot (v. 1.0.0)`
  * `lemon (v. 0.4.3)`
* The analysis supposes the use of 16 processor cores.

#### Running analysis
Before running the analysis be sure to generate the mothur compatible SILVA reference file and indicate in the Makefile its location. The manuscript and supplementary information can be regenerated on a Linux computer by running the following commands:
```
git clone https://github.com/MicrobesRovinj/Korlevic_SeawaterDynamics_x_2021.git
cd Korlevic_SeawaterDynamics_x_2021/
make submission/manuscript.pdf
```
If something goes wrong and the analysis needs to be restarted run the following command from the project home directory before rerunning the analysis:
```
make clean
```
