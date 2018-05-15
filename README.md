# pdb2entropy

PDB2ENTROPY computes entropies from conformational ensembles in PDB format.
A utility program (seq2unfent) is also provided to estimate the entropy of the 
unfolded state of a protein from its sequence.

PDB2ENTROPY program takes in input two files containing:

1) conformational ensembles of the same molecule(s) in PDB format with conformational samples between the MODEL and ENDMDL lines;  
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
PDB2ENTROPY and PDB2TRENT: conformational and translational-rotational entropy from molecular ensembles  
J. Chem. Inf. Model. (submitted)

============================================================================

In the built-in superposition tool, routines are suitable modifications of
those written by D.L. Theobald, therefore, if you use these routines, please
cite:  
Theobald, D. L.  
Rapid calculation of rmsds using a quaternion-based characteristic polynomial.   
Acta. Crystallogr. A, 61, 478–480, 2005.  

============================================================================



COMPILATION:

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
-n k (list entropies based on the first k neighbours (20 default))   
-ne k (entropy estimation based on the kth neighour (10 default))   
-mi (compute entropy using MIST)   
-kmi k (torsions are grouped within the same residue in groups of up to k neighbours. Mutual information among groups will involve at most 2k torsions. (default k 1))   
-s k (use only one snapshot every k snapshots)  
-c X (cutoff distance (Angstrom) for considering mutual information between two groups (default 6.0 A))   
-mr X (minimum resolution (radians) assumed to avoid log(0), 5e-4 default)   
-nt X (number of threads to be used, if less than 1 the program finds the number of threads available)   
-nort (do not superpose all structures to the first one)   
-wp pdb_file (write superimposed strcutures in pdb_file)   
-l (list computed angles)   
-v (verbose mode)   

Usage examples:

--- compute entropy for each residue and sum all entropies. Use 8 threads for parallel computation, do not superpose all structures on first one:

./pdb2entropy sample.pdb tors_next_def.dat sample.out -nt 8 -nort

--- compute entropy using MIST. Consider only first order mutual information for torsions closer in space than 8.0 A. Use 8 threads for parallel computation. Superpose all structures on the first one:

./pdb2entropy sample.pdb tors_next_def.dat sample_mi_1.out -c 8.0 -mi -kmi 1  -nt 8   

--- compute entropy using MIST. Consider mutual information among groups of up to 2 torsion angles, when at least a couple of torsions is closer in space than 8.0 A. Superpose all structure on the first one and write rotated-translated structures in sample_sup.pdb. Use all available threads: 

./pdb2entropy sample.pdb tors_next_def.dat sample_mi_2.out -c 8.0 -mi -kmi 2  -nt 0 -wp sample_sup.pdb 

OUTPUT

For both maximum information spanning tree (MIST) and residue computations, entropies based on the kth (10th by default) nearest neighbour are listed in the output file. First the total entropy and then the entropies for each residue are reported. For MIST calculations the mutual information is divided equally between the residues, when it involves more than one residue. 

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

and finally the sum of all groups entropies minus MIST mutual informations:
- the k^th nearest neighbour (k = 1..20 by default)
- the total entropy value (in R units)

NOTE
The list of entropies using different k neighbours highlights how good is sampling. The entropy should be fairly independent of k, except perhaps for the first neighbours, which have higher variance. 
We found effective, with some thousands conformational samples for proteins to use the MIST approach with m = 1 or 2.
The total entropy corresponding to k = 10 provides typically a stable entropy estimate.  
sample.pdb is provided here only for demonstrative purposes, and to reduce the computational time. Many more conformational samples (in the range of thousands) are needed for accurate estimations of entropy.   

More example files are available in the Download menu at biophysics.uniud.it  

Utility program SEQ2UNFENT

The program seq2unfent takes in input a sequence in fasta format and outputs
the single residues entropies and total entropy of the unfolded state based
on the table (Tab. 1) of Fogolari et al., PLOS One, 10, e0132356, 2015.

The utility program seq2unfent is compiled with:  

cc seq2unfent.c -o seq2unfent -lm   

It is run issuing the command:   

./seq2unfent file_sequence.fas file_entropy.out   

where file_sequence.fas is the file containing the sequence in fasta format and file_entropy.out is the output file containing single residues and total entropies of the unfolded state.

Usage example:  

./seq2unfent b2m.fas b2m_unfolded_entropy.out  
