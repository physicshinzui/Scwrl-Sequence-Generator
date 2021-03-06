# Sequence Generator For Scwrl4
(Scwrl)[http://dunbrack.fccc.edu/SCWRL3.php/#availability] is a software that enables us to mutate amino-acid sequences. 
However, Scwrl does not take into account solvated water.
This script gives a sequence file that is fed to Scwrl, and the sequence file takes into account amino acids that interact with water or other solvent molecules. E.g., lowercase-one-letter amino acids indicate those interacting with solvent, and therefore are fixed in Scwrl rotamer optimisation. 

Quotes:
> Lower-case letters in the sequence indicate that the Cartesian coordinates for the corresponding residues are to be left untouched, and will be treated as steric boundaries only for the other side chains.
http://dunbrack.fccc.edu/SCWRL3.php/#availability

> 1) a new backbone-dependent rotamer library based on kernel density estimates; 
> 2) averaging over samples of conformations about the positions in the rotamer library; 
> 3) a fast anisotropic hydrogen bonding function; 
> 4) a short-range, soft van der Waals atom-atom interaction potential; 
> 5) fast collision detection using k-discrete oriented polytopes; 
> 6) a tree decomposition algorithm to solve the combinatorial problem; 
> 7) optimization of all parameters by determining the interaction graph within the crystal environment using symmetry operators of the crystallographic space group.
http://dunbrack.fccc.edu/SCWRL3.php/#availability

## Requirements
- MDAnalysis >= 2.0  (https://github.com/MDAnalysis/mdanalysis)
- Python >= 3.8.5

## Usage
```
usage: generate_scwrl_seq.py [-h] -i PDB -s FASTA -mr MUTATING_RESIDUE -m MUTATION

optional arguments:
  -h, --help            show this help message and exit
  -i PDB, --pdb PDB     PDB file whose one amino-acid residue mutates.
  -s FASTA, --fasta FASTA
                        fasta file corresponding to the PDB given by `-i`
  -mr MUTATING_RESIDUE, --mutating_residue MUTATING_RESIDUE
                        Usage: 3letter + chainID + residue number in PDB. E.g., ASNE501
  -m MUTATION, --mutation MUTATION
                        One-letter amino acid. E.g., A, K, ...
  -o OUTPUT, --output OUTPUT
```

## Example 
`./generate_scwrl_seq.py -i test/6m0j.pdb -s test/6m0j.fasta -mr ASNE501 -m Y -o scwrl_seq.inp`