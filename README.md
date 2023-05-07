# INTRODUCTION

FASPR is a fast and accurate program for protein side-chain packing, which is an important step in conventional, energy-based protein structure prediction and protein design.

# USAGE

<b> $path/FASPR -i input.pdb -o output.pdb [-s sequence.txt] </b>

-i: the input backbone in pdb format for packing. There should be no missing main-chain atoms (N, CA, C and O). If side-chain atoms are included in input.pdb, they are ignored by FASPR.

-o: the repacked model output in pdb format by FASPR. The residue positions are kept identical to the input pdb file.

-s (optional): an input sequence to be repacked on the input protein backbone. The sequence should be written as one-single line of one-letter alphabet of amino acid types, e.g., ACDEFGHIKLMNPQRSTVWYYWVTSRQPNMLKIHGFEDCA. Only 20 canonical amino-acid types are allowed in the sequence. When the input sequence is the same as the one extracted from the pdb structure, the amino acid side-chain conformations are repacked only. When the input sequence is different, mutations will be introduced. Thus, <b>FASPR can be used to build mutant models efficiently</b>. If you want to fix the conformation of some residues during packing, you can specify them in lower-case letters while keeping other residues in upper-case, e.g., acdefghiklmnpqrstvwyYWVTSRQPNMLKIHGFEDCA.

Please remember to put the executable FASPR program and the binary Dunbrack rotamer library 'dun2010bbdep.bin' in the same directory. Otherwise, the program will report "error! cannot find rotamer library dun2010bbdep.bin". Also, do not change the name of 'dun2010bbdep.bin' because it is hard-coded in the source code.

# INSTALLATION
We recommend users to download the FASPR source-code package to your computer and build the FASPR executable on your own. After downloading and unzipping the package, change into the $path/FASPR/ directory and run "<b>g++ -O3 --fast-math -o FASPR src/*.cpp</b>" if you are working on UNIX or Linux. For Mac users, use "-fast-math" or ignore it. If you are working on the Windows system, you need to install the g++ compiler first.

# COPYRIGHT & CONTACT
Copyright (c) Xiaoqiang Huang. FASPR is free to academic users. For suggestions, please contact xiaoqiah@umich.edu or xiaoqiah@outlook.com.

# REFERENCES
1. Huang x, Pearce R, Zhang Y, FASPR: an open-source tool for fast and accurate protein side-chain packing. Bioinformatics (2020) 36: 3758-3765.
2. Huang X, Pearce R, Zhang Y. Toward the Accuracy and Speed of Protein Side-Chain Packing: A Systematic Study on Rotamer Libraries. J. Chem. Inf. Model. 2020; 60:410-420.
