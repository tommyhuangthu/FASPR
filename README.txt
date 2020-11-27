##############################################################################

FASPR is a method for fast and accurate protein side-chain packing, which is 
an important step in protein structure prediction and protein design.

Copyright (c) Xiaoqiang Huang
Yang Zhang Lab
Dept. of Computational Medicine and Bioinformatics
Medical School
University of Michigan
Email: tommyhuangthu@foxmail.com, xiaoqiah@umich.edu

##############################################################################
                                  FAQ
1. How to use the program?
Usage: ./FASPR -i input.pdb -o output.pdb [-s sequence.txt]
-i: input a pdb backbone for packing. There should be no missing main-chain 
    atoms (N, CA, C and O). If side-chain atoms are included in input.pdb,
    they are ignored by FASPR.
-o: the repacked structure output by FASPR (in standard PDB format). The 
    residue positions are kept identical to the input pdb file
-s: (optionally) pack a new sequence on the backbone of the input.pdb file. 
    The sequence should be written as one-single line of one-letter alphabet 
    of amino acid types, e.g.:
    ACDEFGHIKLMNPQRSTVWYYWVTSRQPNMLKIHGFEDCA
    for the input sequence, only 20 canonical amino acids are allowed.

    When the input sequence is the same as the one extracted from the PDB 
    structure, the amino acid side-chain conformations is repacked. When the 
    input sequence is different from that extracted from PDB, mutations will 
    be introduced. Therefore, FASPR can be used to construct mutant models 
    very efficiently.

    If you want to fix the conformation of some residues during packing, you
    can specify those residue in low-case letters and the others in upper-case 
    letters, e.g.:
    acdefghiklmnpqrstvwyYWVTSRQPNMLKIHGFEDCA

#### NOTE ####
Please remember to put the executable FASPR program and the binary Dunbrack 
rotamer library 'dun2010bbdep.bin' in the same directory. Otherwise, the program 
will prompt "error! cannot find rotamer library dun2010bbdep.bin". Do not change
the name of 'dun2010bbdep.bin' because it is hard coded in the source code.

2. How to build the program for speed?
Build the executable program from the source code using the following command:
g++ --ffast-math -O3 -o FASPR src/*.cpp
or (if the above command does not work):
g++ -O3 -o FASPR src/*.cpp

3. In which OS can I run FASPR?
FASPR can be run in any OS if it is buit following the way as mentioned above.

4. Please report bugs to Dr. Xiaoqiang Huang at the above emails.

5. FASPR is open-source software under the MIT license.

##############################################################################
                                REFERENCE
If FASPR is useful to your work, please cite:

[1] Xiaoqiang Huang, Robin Pearce, Yang Zhang, FASPR: an open-source tool for fast 
and accurate protein side-chain packing. Bioinformatics (2020) 36: 3758-3765.
[2] Xiaoqiang Huang, Robin Pearce, Yang Zhang. Toward the Accuracy and Speed of 
Protein Side-Chain Packing: A Systematic Study on Rotamer Libraries. 
J. Chem. Inf. Model. 2020; 60:410-420.
