# Overview
This python routine is able to calculate parameters associated with a two-state folding model with data obtained from a measured observation such as absorbance or 
fluorescence.


# Background 
Folding curves for folding and unfolding using a two state model. For more information on the two-state model, read the following papers:

 1. "Unfolding free energy changes determined by the linear extrapolation method. 1. Unfolding of phenylmethanesulfonyl .alpha.-chymotrypsin using different denaturants." Marcelo M. Santoro , D. W. Bolen, Biochemistry, 1988, 27 (21), pp 8063–8068, DOI: 10.1021/bi00421a014
 2. "Folding of horse cytochrome c in the reduced state." Bhuyan AK, Udgaonkar JB, J Mol Biol. 2001 Oct 5;312(5):1135-60.

## Equations 
There is a discrepancy between several publications for the folding equetion.

From Bhuyan's paper:

Sobs=(Cf +mf[D]+Cu+mu[D]exp((-Delta G +mg[D])/(RT))/(1+exp((-DeltaG +mg[D])/(RT)))

From Santoro's paper:

DeltaEpsilon=[(Cf+mf[D])+(Cu+mu[D])*exp-(Delta G/RT +mg[D]/RT)]/[1+exp-(Delta G/RT +mg[D]/RT)]


These two equations are inconsistent with one another. 

For the fitting routine, the following equation is used:

Sobs=(Cf +mf[D]+(Cu+mu[D])exp((-Delta G +mg[D])/(RT))/(1+exp((-DeltaG +mg[D])/(RT)))

Using this equation is able to reproduce the numbers found in the papers. 


# Running the program
For proof of principle, the program can be run to show a typical output with the command using the default data.

`$ python folding-curves.py`

This will show a plot and also ask if you want to save the parameters to CSV files.

To run the routine with a specified CSV file

`$ python folding-curves.py my-file.csv`

