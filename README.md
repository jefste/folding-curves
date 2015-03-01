# folding-curves
Folding curves for folding and unfolding using a two state model
ssdfsdfsdfsdfdf
## Equations 
There is a discrepancy between several publications for the folding equation.

From 'Folding of Reduced Horse Cytochrome c' Bhuyan and Udgaonkar JMB, doi:10.1006/jmbi.2001.4993
Sobs=(Cf +mf[D]+Cu+mu[D]exp((-Delta G +mg[D])/(RT))/(1+exp((-DeltaG +mg[D])/(RT)))

From 'Unfolding Free Energy Changes Determined by the Linear Extrapolation method...', Santoro and Bolen, Biochemistry, 
doi: 10.1021/bi00421a014

DeltaEpsilon=[(Cf+mf[D])+(Cu+mu[D])*exp-(Delta G/RT +mg[D]/RT)]/[1+exp-(Delta G/RT +mg[D]/RT)]

For the fitting routine, the following equation is used:
Sobs=(Cf +mf[D]+(Cu+mu[D])exp((-Delta G +mg[D])/(RT))/(1+exp((-DeltaG +mg[D])/(RT)))

