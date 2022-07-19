# GFIM
Generate Full Interaction Model

This code is design to : 
- generate a quadratic design matrix with simple, interaction and pure quadratic effects by default the simple effect matrix is binomial.
- generate a sparse parameter model by breaking down the non-zero coordinates into simple, interaction and quadratic effects

There is no argument to give, but you can modify the code to adjust : 
- Number of SNPs
- Number of individuals
- h2
- ratio 0/1 (binom rule)
- number of signals
