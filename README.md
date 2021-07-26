# Nonadditivity analysis


## Synposis


A program to find key complex patterns in SAR data


--------------------


## Installation

The program has been tested on Python 3.6.

You will need a copy of the RDKit cheminformatics toolkit, available
from http://rdkit.org/. The easiest way is to install via PyPI with
`pip install rdkit-pypi`.

Install directly from source with:

```bash
$ pip install git+https://github.com/KramerChristian/NonadditivityAnalysis.git
```

Install the code in development mode with:

```bash
$ git clone git+https://github.com/KramerChristian/NonadditivityAnalysis.git
$ cd NonadditivityAnalysis
$ pip install -e .
```

The path to mmpdb has to be set on line 44 of the Nonadditivity analysis 
code. If a special salt clean-up is required, the path to the salt definitions
can be set on line 43.

## How to run the program and get help

The code runs as a simple command-line tool. Command line options are printed via

```shell
$ python -m nonadditivity -h
```

## Example usage

Using the test files supplied, an example run can be

```shell
$ python -m nonadditivity -in hERG_ChEMBL.txt -delimiter tab -series_column ASSAY_CHEMBLID -props PCHEMBL_VALUE -units nM
```

#### Input file format

IDENTIFIER [sep] SMILES	[sep] DATA
...

where [sep] is the separator and can be chosen from tab, space, comma, and 
semicolon.


------------------


## Publication

If you use this code for a publication, please cite
Kramer, C. Nonadditivity Analysis. J. Chem. Inf. Model. 2019, 59, 9, 4034–4042.

https://pubs.acs.org/doi/10.1021/acs.jcim.9b00631


-----------------


## Background

The overall process is:

  1) Parse input:
     - read structures
     - clean and transform activity data
     - remove Salts

  2.) Compute MMPs

  3.) Find double-transformation cycles

  4.) Write to output & calculate statistics


#### 1) Parse input

Ideally, the compounds are already standardized when input into nonadditivity 
analysis. The code will not correct tautomers and charge state, but it will 
attempt to desalt the input.

Since Nonadditivity analysis only makes sense on normally distributed data, the
input activity data can be transformed depending on the input units. You can choose
from "M", "mM", "uM", "nM", "pM", and "noconv". The 'xM' units will be transformed
to pActivity wiht the corresponding factors. 'noconv' keeps the input as is and does
not do any transformation.

For mulitplicate structures, only the first occurence will be kept.


#### 2) Compute MMPs

Matched Pairs will be computed based on the cleaned structures. This is done by a
subprocess call to the external mmpdb program. Per default, 20 parallel jobs are used
for the fragmentation. This can be changed on line 681.


#### 3) Find double-transformation cycles

This is the heart of the Nonadditivity algorithm. Here, sets of four compounds that are
linked by two transformations are identified. For more details about the interpretation
see publication above.


#### 4) Write to output and calculate statistics

Information about the compounds making up the cycles and the distribution of 
nonadditivity is written to output files. [...] denotes the input file name.
The file named 

"Additivity_diffs"[...]".txt"

contains information about the cycles and the Probability distribution


The file named

"Additivity_diffs"[...]"_perCompound.txt"

contains information about the Nonadditivity aggregated per Compound across all cycles
where a given compound occurs. 


The file named

"Additivity_diffs"[...]_c2c.txt

links the two files above and can be used for examnple for visualizations in SpotFire.


--------------------


## Copyright

The NonadditivityAnalysis code is copyright 2015-2019 by F. Hoffmann-La
Roche Ltd and distributed under the 3-clause BSD license (see LICENSE.txt).

