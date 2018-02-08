# pdb2entropy
PDB2ENTROPY: entropy calculation from conformational ensembles

PDB2ENTROPY computes entropies from conformational ensembles in PDB format.

The program takes in input two files containing:

1) conformational ensembles of the same molecule(s) in PDB format;
2) definitions of torsion angles (a default file is provided with the program
where additional user's definitions can be easily implemented).

The program outputs residue-based entropies or maximum-information spanning
tree entropies, together with information related to the specific method used.
