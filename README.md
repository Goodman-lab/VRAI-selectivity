
# VRAI-selectivity REAMME.md

This document contains the manual for VRAI-selectivity.py and VRAI-multi.py (extension to the VRAI-selectivity.py script (version > v7))


===============================================================

CONTENTS

VRAI-selectivity Manual
1) Requirements and Setup
2) Usage
3) Example Usage

VRAI-multi Manual
1) Requirements and Setup
2) Programme Documentations
3) Examples

===============================================================

# References for the underlying theory

1. S. Lee and J. M. Goodman,VRAI-selectivity: calculation of selectivity beyond transition state theory, Org. Biomol. Chem., 2021, 19, 3940–3947. (https://doi.org/10.1039/D1OB00234A)
2. 	C. C. Lam, J. M. Goodman, Reaction dynamics as the missing puzzle: the origin of selectivity in oxazaborolidinium ion-catalysed reactions, 	Chem. Sci., 2023,14, 12355-12365 (https://doi.org/10.1039/D3SC03009A)
   
      (Feel Free to drop Ching Ching an email if there is any question on the work from this paper: ccl59@cam.ac.uk)

# VRAI-selectivity Manual: 
===============================================================

VRAI-selectivity Manual

Version 1.0

Copyright (c) 2020 Sanha Lee, Jonathan M. Goodman
University of Cambridge

This documents is adapted from 'VRAI-Selectivity: better 
calculation of selectivity than transition state theory' 
Supporting Information

===============================================================

CONTENTS
1) Requirements and Setup
2) Usage
3) Example Usage

===============================================================

1) REQUIREMENTS AND SETUP

The script is currently set up to run on Python version 2.7. The following Python modules must be installed before the script can be run:
- numpy
- argparse
- os
- math
- sys
- rdkit
- bisect
- pprint

The python script is designed to analyse Gaussian output files from Gaussian 16 only.

The script will only recognise the Gaussian output file if the frequency calculation is run independently. Do not run the calculation with 'opt' and 'freq'  keywords at the same time. The following frequency keyword must be included in the calculation 'freq=(hpmodes,intmodes)'.  

The python script requires the Gaussian output files to be in .mol file format.
It is possible to convert the frequency output files to .mol file using Gaussview via File -> Convert Files.


2) CORRECT USAGE

The .mol files and .out files must be inputed in the following order.

python VRAI-selectivity.py TS1molfile.mol INTmolfile.mol P1molfile.mol P2molfile.mol TS1freq.out TS2freq.out [OPTIONAL TS2Afreq.out] [OPTIONAL TS2Bfreq.out] [OPTIONAL -I]

TS1molfile.mol - geometry of the first transition state in .mol file format
INTmolfile.mol - geometry of the intermediate in .mol file format 
P1molfile.mol - geometry of the first product in .mol file format
P2molfile.mol - geometry of the second product in .mol file format
TS1freq.out - Gaussian16 frequency calculation output file for the first transition state
INTfreq.out - Gaussian16 frequency calculation output file for the intermediate
TS2Afreq.out - Gaussian16 frequency calculation output file for TS2A (optional)
TS2Bfreq.out - Gaussian16 frequency calculation output file for TS2B (optional)
-I - optional argument for intermediate activation
-S - optional argument for TST calculation

The analysis of reactions proceeding on the potential energy surface with an intermediate should use '-I' for intermediate activation.
- When -I option is used, the TS2A and TS2B Gaussian frequency output files must be provided (the algorithm will return an error message if these are not provided)
- When -I option is not used, do not provide TS2A or TS2B Gaussian frequency output files (the algorithm will return an error message if these are provided)
- The TS2A must correspond to the P1 forming TS, connecting INT with P1
- The TS2B must correspond to the P2 forming TS, connecting INT with P2
The algorithm will read in the TS2 frequency output files and extract the free energy of the molecule. The free energy difference between TS1 and TS2 will be used to decide whether the algorithm uses TST or dynamic analysis for the selectivity calculation.

When '-S' activation is used, the algorithm will simply output the product distribution calculated from TST. This use only requires the TS2A and TS2B Gaussian 16 output files as the input. The correct execution of the code on the command line would be:

python VRAI-selectivity.py TS2Afreq.out TS2Bfreq.out -S

For bifurcation analysis, the INTmolfile should be replaced by TS2molfile which is the geometry of the second transition state. Do not use intermediate activation.

The key output will be in the following format:
- files corresponding to the major and minor product are stated
- the two μ and λ values which decides the alignment type are given 
- the length of vector |g| and the angle φ (angle between the reaction vector and vector g_ is stated

The algorithm will create a log file and print out
- whether the rdkit object was successfully created
- vectors p1_, p2_ and the imaginary eigenvector
- φ: angle between the imaginary eigenvector, a_, and vector g_ (displacement of TS2 from TS1)
- dot products, angles and magnitude of the key vectors
- uncommon bonds the algorithm has used for the selectivity prediction
- returns a warning sign for large |g_| and when the predicted selectivity is low


3) EXAMPLE USAGE

The following files are provided for this section
- testP1freq.mol
- testP2freq.mol
- testTS1freq.mol
- testTS1freq.out
- testINTfreq.mol
- testINTfreq.out
- testTS2A.out
- testTS2B.out

Test the code on the example Gaussian output files given: 
python VRAI-selectivity_v2.py testTS1freq.mol testINTfreq.mol testP1freq.mol testP2freq.mol testTS1freq.out testINTfreq.out testTS2Afreq.out testTS2Bfreq.out -I

The Gaussian output files were generated by optimising the reaction published by Singleton et al in 2006 (Ussing, B. R.; Hang, C.; Singleton, D. A., J. Am. Chem. Soc. 2006, 128, 7594–7607). The correct execution should print out the following output.

**** Analysis Completed ****
Major product is testP2freq.mol
Minor product is testP1freq.mol

mu1_ = -0.0441868111245
mu2_ = 0.0484860277906
lambda1_ = 0.432580161859
lambda2_ = 0.438546804793
|g_| = 0.483733089154
phi = 6.62625149713

The algorith will now proceed to estimate the major and minor product ratios

Product Ratio Calculation Completed:
Major Product : Minor Product ratio
82.9 : 30.1

****************************

The algorithm should also create a log file called 'VRAI-selectivity_testTS1freq.log' with all the key information.



# VRAI-multi Manual:

VRAI-multi.py is an extension to the VRAI-selectivity.py script (version > v7). VRAI-multi automate the process of VRAI-selectivity analyses for treating systems with complex PES.


===============================================================

CONTENTS
1) Requirements and Setup 
2) Programme Documentations
3) Examples

===============================================================


## 1. Requirements and Setup

### Inputs 
The input to VRAI-multi is the path to the folder that contains the required files for the VRAI selectivity analyses. For the script to correctly identify and process the files, the following file naming formats must be adopted: 

-	Geometry of the first TS in mol file must have a ‘TS1.mol’ suffix. 
-	Gaussian16 frequency calculation output file for the first TS must adopt the same filename as the corresponding mol file but with a ‘.out’ suffix.
```
#Example: 
ol17ryu_e_int2_SSR_11_TS1.mol and ol17ryu_e_int2_SSR_11_TS1.out
```
-	Geometry of the intermediate in .mol file must have a ‘int.mol’ suffix. 
-	Gaussian16 frequency calculation output file for the intermediate must adopt the same filename as the corresponding mol file but with a ‘.out’ suffix. 
```
#Example: 
ol17ryu_e_int2_SSR_6_int.mol and ol17ryu_e_int2_SSR_6_int.out
```
-	Gaussian16 frequency calculation output files for the second TS must have a ‘TS2.out’ suffix. (optional)
-	Geometry of the product in .mol file must adopt the same filename as the corresponding TS2 with an additional ‘prod.mol’ as the suffix.
```
#Example: there are three different possible products via pathways that share a common intermediate 
ol17ryu_e_int2_SSR_C_1_TS2_prod.mol
ol17ryu_e_int2_SSR_C_1_TS2.out
ol17ryu_e_int2_SSR_H_2_TS2_prod.mol
ol17ryu_e_int2_SSR_H_2_TS2.out
ol17ryu_e_int2_SSR_O_1_TS2_prod.mol
ol17ryu_e_int2_SSR_O_1_TS2.out
```

- Single point energy calculation files must adopt the same filename as the corresponding TS1, int and TS2 structures but with a '_spe.out’ suffix. (optional)
```
#Example: 
ol17ryu_e_int2_SSR_O_1_TS2_spe.out for ol17ryu_e_int2_SSR_O_1_TS2.out
```
The designated folder must contain a minimum of one set of TS1 files, one set of INT files and two product mol files. It is acceptable to include additional sets of product or TS1 files, but the number of intermediate files in the folder should not exceed one set (ie a .mol and .out file).

### Environment

The scripts are currently set up to run on Python version 3.8. The following Python modules must be installed before the script can be run:

- numpy
- pandas 
- rdkit

The follow packages are included in the Python Standard Library: 
- argparse
- os
- math
- sys
- bisect
- pprint
- itertools
- datetime 

The python script is designed to analyse Gaussian output files from Gaussian 16 only.

The script will only recognise the Gaussian output file if the frequency calculation is run independently. Do not run the calculation with 'opt' and 'freq'  keywords at the same time. The following frequency keyword must be included in the calculation 'freq=(hpmodes,intmodes)'.  

The python script requires the Gaussian output files to be in .mol file format.
It is possible to convert the frequency output files to .mol file using Gaussview via File -> Convert Files.


## 2. Programme Documentations


### class VRAI_multi

### VRAI_multi.__init__(self, file_info) 


file_info: string; the path to the folder that contains all the input files or the path to the csv that contains the name of the input files 


### VRAI_multi.autopipe(self,path='', intermediate_option=True, weight_option=True, TST_option=False, getcsv=False, with_spe=False, temperature = 298.0, energy_cut_off=5 )

#### Parameters 

path: string; the path to the folder that contains all the input files if the path to a csv file is given for file_info

intermediate_option: bool; default = True; the argument for intermediate activation - should set to True if TS2s are presented on the PES 

weight_option:  bool; default = True; Alignment weight option - True = Alignment weight option is set to atomic weights; False = Alignment weight option is set equal for all atoms

TST_option: bool; default = False; apply TST theory when TS2 energy is higher than -(energy_cut_off) kJ/mol from TS1

energy_cut_off: float; default = 5; energy cut-off for activating product ratio calculation based on TST; unit: kJ/mol

with_spe: bool; default = False; take account of single point energy outputs at a higher level of theory when calculating the free energy

getcsv: bool; default = False; return the output as .csv file

temperature: float; default = 298.0


#### Attributes


raw_result_df: pd.dataframe; raw results from VRAI-selectivity.py

processed_result_df: pd.dataframe; The output of the pipeline is a data frame that contains the product percentages, which are calculated by using different products as the reference in the calculation. It is crucial to verify the consistency of the product percentages obtained from different sets of results to accept them as reliable outcomes.


## 3. Examples

```python
# import the script
from VRAI_multi import VRAI_multi

```


### Example 1: 
- With file_info = the path to the folder that contains all the input files

```python

test=VRAI_multi('./test_file2/')

# with_spe=True: take account of single point energy outputs at a higher level of theory when calculating the free energy
test.autopipe(with_spe=True)

Raw_result = test.file_df
Result = test.processed_result_df

```

### Example 2:
- With file_info = the path to the csv that contains the name of the input files


```python

test2=VRAI_multi('test_file.csv')
test2.autopipe(path='./test_file2/',with_spe=True)


```

