# SOLkq
Lanthanide Stark energy levels calculator

This program calculates the Bkq parameters of the Ligand Field Hamiltonian of 4f^n Lanthanide ions (Ce3+ - Yb3+) with optically active 4f electrons,
according to the Simple Overlap model developed by Oscar Malta (O.L. Malta, Chem. Phys. Lett. v.87 (1), 1982)

The input file should be an .xyz geometry file (following the mol2 style, first line is the number of atoms, second line additional information, and geometry starts on third line - important!!!).

The file should contain only the first coordination sphere and the Ln3+ ion should be the first atom - important!!!. In case of doubt, check the provided examples of .xyz files (one with integer J, Eu3+, one with half-integer J, Er3+).

You need to provide the charges of the atoms (Write them separated by spaces, such as for instance, 4 atoms: 1 0.5 0.2 1.5 1). I strongly suggest obtaining them through the Omega(lambda) parameters fit. One can do that through the JOYSpectra web platform. For the example tests, use just 1 for all charges to get the hang of the program

This program does not need installing. Simply copy the .exe file to the folder of your .xyz file. The output will be printed in the same folder

Supported ligands: N, O, F, P (Not for Ho3+), S, Cl, Se (not for Ho3+, Tm3+, and Lu3+), Br (not for Er3+), I
Supported Lanthanides (III): Ce3+, Pr3+, Nd3+, Pm3+, Sm3+, Eu3+, Gd3+, Tb3+, Dy3+, Ho3+ (not for P, and Se), Er3+ (not for Br), Tm3+ (not for Se), and Yb3+

Output: Bkq values and Ligand Field splitting for selected (2S+1)L(J) states.

Note1: Half-Integers J states need to be split into Krammers doublets
Note2: Matrix elements calculated in the RUSSEL-SAUNDERS coupling scheme, with LS-Coupling matrix elements taken from the tables of Nielson and Koster.

If used, please reference this program. The paper should soon follow, future DOI will be added.
PSA: I'm a Chemist and not a Programmer, this software is by no means optimized, but it calculates what it needs to.
