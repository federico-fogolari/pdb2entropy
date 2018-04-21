# pdb2entropy

PDB2ENTROPY computes entropies from conformational ensembles in PDB format.

The program takes in input two files containing:

1) conformational ensembles of the same molecule(s) in PDB format;
2) definitions of torsion angles (a default file is provided with the program
where additional user's definitions can be easily implemented).

The program outputs residue-based entropies or maximum-information spanning
tree entropies, together with information related to the specific method used.

CONTACT:  

Federico Fogolari  
Dipartimento di Scienze Matematiche, Informatiche e Fisiche  
Universita' di Udine  
Via delle Scienze 206  
33100 Udine - Italy  
Tel ++39 0432 494320    
E-mail federico.fogolari@uniud.it  


REFERENCE:

Please cite:  
F. Fogolari, O. Maloku, C.J. Dongmo Foumthuim, A. Corazza, G. Esposito  
PDB2ENTROPY and PDB2TRENT: entropy from conformational ensembles  
J. Chem. Inf. Model. (submitted)

In the built-in superposition tool, routines are suitable modifications of
those written by D.L. Theobald, therefore, if you use these routines, please
cite:  
Theobald, D. L. (2005).   
Rapid calculation of rmsds using a quaternion-based characteristic polynomial.   
Acta. Crystallogr. A, 61, 478–480.  


INSTALLATION:

The program is compiled with: 

- if OpenMP is installed:  
cc pdb2entropy.c -o pdb2entropy -lm -fopenmp

- otherwise:  
cc pdb2entropy.c -o pdb2entropy -lm 

FORMAT OF TORSION-ADJACENCY FILE:

The file tors_next_def.dat is provided as a sample file.   
There are two type of records in this file:

- Adjacency record
NEXT  residue_name  atom_name1  atom_name2  cutoff_distance  

A residue is adjacent to the next if its atom1 is within cutoff_distance from atom_name2 of the next residue (whatever its name).  
The next residue is identified by the program by checking residue number and insertion character within the same segment and chain. 

- Torsion record
TORS residue_name torsion_name atom_name1 atom_name2 atom_name3 atom_name4 symmetry

The torsion_name torsion angle for residue_name is the dihedral defined by the four atom_name in sequential order. Symmetry is an integer defining the rotational symmetry of the torsion angle (e.g. 2 for rotation of a carboxylic group, 3 for rotation of a methyl group).  
The signs + and - indicate atoms on the following and previous adjacent residue, respectively.

- Angle records (not present in the sample tors_next_def.dat file) in the form:
ANG residue_name angle_name atom_name1 atom_name2 atom_name3  
are also recognized.

The sample file tors_next_def.dat contains data for proteins and nucleic acids.

RUNNING PDB2ENTROPY

./pdbentropy without arguments will print options available

./pdb2entropy expects at least three arguments:
 - the name of the input pdb_file 
 - the name of the torsion-adjacency definition file
 - the name of the output file

Other options are listed hereafter

Usage:
./pdb2entropy pdb_infile def_infile outfile [Options]
Options:
-n (max k neighbours for listing entropies (20 default))
-mi (compute entropy from MIST)
-kmi k (compute entropy from MIST considering mutual information among groups of k variables (default k 1))  
-c X (cutoff distance (Angstrom) for MI pair filtering)
-mr X (minimum resolution (in radians) assumed to avoid log(0), 5e-4 default)
-nt X (number of threads to be used, if less than 1, e.g. with -nt 0, the program finds the number of threads available)
-l (lists all defined torsion angles as rows with values for each conformation)
-nort (do not superpose all structures on the first one)
-wp pdb_file (write superimposed structures in pdb_file)
-v (verbose mode)

Usage examples:

--- compute entropy for each residue and sum all entropies. Use 8 threads for parallel computation, do not superpose all structures on first one:

./pdb2entropy sample.pdb tors_next_def.dat sample.out -nt 8 -nort

--- compute entropy using MIST. Consider only first order mutual information for torsions closer in space than 8.0 A. Use 8 threads for parallel computation. Superpose all structure on the first one:

./pdb2entropy sample.pdb tors_next_def.dat sample_mi_1.out -c 8.0 -mi -kmi 1  -nt 8   

--- compute entropy using MIST. Consider mutual information among groups of up to 2 torsion angles, when at least a couple of torsions is closer in space than 8.0 A. Superpose all structure on the first one and write rotated-translated structures in sample_sup.pdb. Use all available threads:
./pdb2entropy sample.pdb tors_next_def.dat sample_mi_2.out -c 8.0 -mi -kmi 2  -nt 0 -wp sample_sup.pdb 

OUTPUT

For both maximum information spanning tree (MIST) and residue computations entropies based on the kth (10th by default) nearest neighbour are listed in the output file. First the total entropy and then the entropies for each residue are reported. For MIST calculations the mutual information is divided equally between the residues, when it involves more than one residue. 

A long version of the same file is created with the same name with the extension .long appended to the file name.  This is described hereafter.

1) The output for residue-based entropy lists for each residue 
- the k^th nearest neighbour (k = 1..20 by default)
- the entropy value (in R units)
- the average distance to the k^th nearest neighbour
- the residue (name, number, insertion, chain, segid if present)

and finally the sum of all residue entropies listing:
- the k^th nearest neighbour (k = 1..20 by default)
- the total entropy value (in R units)

2) The output for MIST entropy lists first 
- the groups of torsions. 

Then the output lists for each group  
- the k^th nearest neighbour (k = 1..20 by default)
- the entropy value (in R units)
- the average distance to the k^th nearest neighbour
- the group of torsions
- the residue (name, number, insertion, chain, segid if present) to which
  the group torsions belongs (a detailed list may be obtained by the option -v)

Then the output lists the mutual information entering the MIST among groups:
- the k^th nearest neighbour (k = 1..20 by default)
- the two groups of torsions
- the mutual information value (in R units)

and finally the sum of all groups entropies minus Minimum Spanning tree 
mutual informations:
- the k^th nearest neighbour (k = 1..20 by default)
- the total entropy value (in R units)

NOTE
The list of entropies using different k neighbours highlights how good is sampling. The entropy should be fairly independent of k, except perhaps for the first neighbours, which have higher variance. 
We found effective, with some thousands conformational samples for proteins to use the MIST approach with m = 1 or 2.
The total entropy corresponding to k = 10 provides typically a stable entropy estimate.  
sample.pdb is provided here only for demonstrative purposes, and to reduce the computational time. Many more conformational samples (in the range of thousands) are needed for accurate estimations of entropy.

