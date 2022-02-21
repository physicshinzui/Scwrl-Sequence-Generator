#!/usr/bin/env python3
from MDAnalysis import Universe
import argparse

def AAresidue2index(univ, verbose=False):
    """
    Args: 
        univ: Universe object of MDAnalysis
        verbose: if true, print altLoc (Default is False)

    Return: 
        resi2indx: a dictionary that maps a string (3-letter AA + chainID + resid; e.g., "ASNE501") to 0-start consecutive index (note: the index is consecutive even for multiple chains.)    

    >>> univ = Universe("test/3aa.pdb")
    >>> AAresidue2index(univ)
    {'ALASYSTEM2': 0, 'ASPSYSTEM3': 1, 'TRPSYSTEM4': 2}
    >>> univ = Universe("test/6m0j.pdb")
    >>> AAresidue2index(univ)['SERA19'], AAresidue2index(univ)['GLYE526'], AAresidue2index(univ)['ASNE501']
    (0, 790, 765)
    """
    protein = univ.select_atoms("protein and name CA")
    
    icou = 0
    resi2indx = {}
    for aa in protein:
        
        if aa.altLoc == '':
            resi2indx[aa.resname+aa.segid+str(aa.resid)] = icou
            icou += 1      

        elif aa.altLoc == 'A':
            resi2indx[aa.resname+aa.segid+str(aa.resid)] = icou
            if verbose: print(f"{icou}, {aa.resname}, {aa.resid}: Only altLoc {aa.altLoc} was considred.")
            icou += 1
    
    return resi2indx

def get_index_of_aa_near_solvent(univ, verbose=False):
    """
    Args:
        univ: Universe object of MDanalysis
    Return: 
        seq_index: Index of amino-acid residues around solvent.
                e.g., Suppose the sequence 'AwD' (w indicates one around solvent). The index of 'w' is 1 (zero-start).
    """
    resi2indx = AAresidue2index(univ)
    near_solvent = univ.select_atoms("(around 4.0 resname HOH ZN* CL* NA*) and protein")
    near_solvent.write('near_solvent.pdb')
    seq_index = []
    for residue in near_solvent.residues:
        res = residue.resname+residue.segid+str(residue.resid)
        seq_index.append(resi2indx[res])
    return seq_index

def oneline_seq(fasta):
    """
    Args:
        fasta: a fasta file

    Return:
        seq: one-liner sequence. Sequences of chains are concatinated.

    >>> oneline_seq("test/3aa.fasta")
    'ADW'
    >>> oneline_seq("test/6m0j.fasta")
    'STIEEQAKTFLDKFNHEAEDLFYQSSLASWNYNTNITEENVQNMNNAGDKWSAFLKEQSTLAQMYPLQEIQNLTVKLQLQALQQNGSSVLSEDKSKRLNTILNTMSTIYSTGKVCNPDNPQECLLLEPGLNEIMANSLDYNERLWAWESWRSEVGKQLRPLYEEYVVLKNEMARANHYEDYGDYWRGDYEVNGVDGYDYSRGQLIEDVEHTFEEIKPLYEHLHAYVRAKLMNAYPSYISPIGCLPAHLLGDMWGRFWTNLYSLTVPFGQKPNIDVTDAMVDQAWDAQRIFKEAEKFFVSVGLPNMTQGFWENSMLTDPGNVQKAVCHPTAWDLGKGDFRILMCTKVTMDDFLTAHHEMGHIQYDMAYAAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKSIGLLSPDFQEDNETEINFLLKQALTIVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEPVPHDETYCDPASLFHVSNDYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLFNMLRLGKSEPWTLALENVVGAKNMNVRPLLNYFEPLFTWLKDQNKNSFVGWSTDWSPYADTNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCG'
    """
    fin = open(fasta, 'r')
    seq=""
    for line in fin:
        if line[0] == ">": continue
        seq += line.rstrip()
    return seq

def substitute_single_aa(seq, aa2indx, residue, mt_single_letter):
    """
    Args:
        seq : a one-line sequence which is returned by the functino `oneline_seq`
        aa2indx: a dict that maps "3 letter AA+chainID+resid" to an 0-start concecutive index. This is returned by `AAresidue2index`.
        residue: a string "3 letter AA+chainID+resid" of an amino acid to mutate.
        mt_single_letter: one-letter an amino acid string 

    Return: 
        sub_seq: a sequence with a single mutation.

    >>> substitute_single_aa('AIK',{'ALAA10':0, 'ILEA11':1, 'LYSA12':2}, 'ILEA11', 'R')
    'ARK'
    >>> u1    = Universe("test/3aa.pdb")
    >>> aa2indx = AAresidue2index(u1)
    >>> seq1     = oneline_seq("test/3aa.fasta")
    >>> substitute_single_aa(seq1, aa2indx, "ASPSYSTEM3", mt_single_letter="K")
    'AKW'
    >>> u2    = Universe("test/6m0j.pdb")
    >>> aa2indx = AAresidue2index(u2)
    >>> seq2     = oneline_seq("test/6m0j.fasta")
    >>> substitute_single_aa(seq2, aa2indx, "ASNE501", mt_single_letter="Y")[aa2indx["ASNE501"]]
    'Y'
    """
    assert len(seq) == len(aa2indx), f"The number of amino acids for seq ({len(seq) }) and pdb ({len(aa2indx)}) are different!, which must be THE SAME."
    sub_seq = ""
    for i, aa in enumerate(seq):
        if i == aa2indx[residue]:
            sub_seq += mt_single_letter.upper()
            #print(" >> "+mt_single_letter.upper()+" << ", end='')
        else:
            sub_seq += aa.upper()
            #print(aa.upper(), end='')

    return sub_seq

def scwrl_seq(seq, ilocs): 
    """
    Args:
        seq: one-line sequence. 
        ilocs: location index (starts from 0) of amino acids to be untatched in the scwrl computation.

    Return: 
        seq_for_scwrl: a sequence to be fed to Scwrl.    
    """
    seq_for_scwrl = ""
    for i, aa in enumerate(seq): 
        if i in ilocs:
            seq_for_scwrl += aa.lower()

        else:
            seq_for_scwrl += aa.upper()

    return seq_for_scwrl

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-i' ,'--pdb', required=True, help='PDB file whose one amino-acid residue mutates.')
    p.add_argument('-s' ,'--fasta', required=True, help = 'fasta file corresponding to the PDB given by `-i`')
    p.add_argument('-mr','--mutating_residue', required=True, help='Usage: 3letter + chainID + residue number in PDB. E.g., ASNE501')
    p.add_argument('-m' , '--mutation', required=True, help='One-letter amino acid. E.g., A, K, ...')
    p.add_argument('-o', '--output', default="scwrl_seq.inp")
    args      = p.parse_args()
    pdb       = args.pdb
    fastafile = args.fasta
    mres      = args.mutating_residue
    mt        = args.mutation
    outseq    = args.output

    u                 = Universe(pdb)
    aa2indx           = AAresidue2index(u, verbose=True)
    aa_near_solv_indx = get_index_of_aa_near_solvent(u)
    seq               = oneline_seq(fastafile)

    print("Input Seq:\n", seq)
    print()

    single_mt_seq = substitute_single_aa(seq, aa2indx, residue=mres, mt_single_letter=mt)
    print("Single-point Mutation:\n", single_mt_seq)
    print()
    scwrl_input_seq = scwrl_seq(single_mt_seq, aa_near_solv_indx)    
    with open(outseq, "w") as fout:
        fout.write(scwrl_input_seq)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    main()
