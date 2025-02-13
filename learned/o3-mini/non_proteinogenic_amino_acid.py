"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
#!/usr/bin/env python3
"""
Classifies: Non‐proteinogenic amino acids
Defined as: “Any amino‐acid that is not naturally encoded in the genetic code of any organism.”
Uses pattern matching for an amino acid backbone plus a comparison with 20 canonical proteinogenic amino acids.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non‐proteinogenic amino acid.
    A molecule is first required to have an amino acid backbone. Then, if it exactly
    matches one of the 20 proteinogenic amino acids (by canonical SMILES), it is not considered non‐proteinogenic.
    Otherwise, if an amino acid backbone is detected and no match with standard amino acids is found,
    it is classified as non‐proteinogenic.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a non‐proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for an amino acid backbone.
    # Pattern 1: Open‐chain alpha‐amino acid: an sp3 carbon with an attached amino group and a carboxylic acid.
    aa_pattern = Chem.MolFromSmarts("[C;H](N)C(=O)[O;H,-]")
    # Pattern 2: Cyclic amino acid (e.g. proline) – note that the amino group is part of a ring.
    proline_pattern = Chem.MolFromSmarts("O=C(O)N1CC[C@H](C1)")  

    # Check if the molecule has an amino acid backbone (either open‐chain or proline type)
    if not (mol.HasSubstructMatch(aa_pattern) or mol.HasSubstructMatch(proline_pattern)):
        return False, "No typical amino acid backbone detected"
    
    # Compute the canonical SMILES for input molecule for standardized comparison
    input_canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    
    # List of representative SMILES for the 20 proteinogenic amino acids.
    # Note: For some amino acids stereochemistry is important so we use one representative.
    proteinogenic_aas = [
        "NCC(=O)O",                              # glycine
        "N[C@@H](C)C(=O)O",                       # L-alanine
        "N[C@@H](C(C)C)C(=O)O",                    # L-valine
        "N[C@@H](CC(C)C)C(=O)O",                   # L-leucine
        "N[C@@H](C[C@H](C)C)C(=O)O",               # L-isoleucine
        "O=C(O)N1CC[C@H](C1)O",                    # L-proline (representative cyclic form)
        "N[C@@H](Cc1ccccc1)C(=O)O",                # L-phenylalanine
        "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",         # L-tryptophan
        "N[C@@H](Cc1ccc(O)cc1)C(=O)O",              # L-tyrosine
        "N[C@@H](CO)C(=O)O",                       # L-serine
        "N[C@@H]([C@H](O)C)C(=O)O",                # L-threonine (one representative)
        "N[C@@H](CS)C(=O)O",                       # L-cysteine
        "N[C@@H](CCSC)C(=O)O",                      # L-methionine
        "N[C@@H](CC(=O)O)C(=O)O",                   # L-aspartic acid
        "N[C@@H](CCC(=O)O)C(=O)O",                  # L-glutamic acid
        "N[C@@H](CC(=O)N)C(=O)O",                   # L-asparagine
        "N[C@@H](CCC(=O)N)C(=O)O",                  # L-glutamine
        "N[C@@H](CCCCN)C(=O)O",                     # L-lysine
        "N[C@@H](CCCNC(=N)N)C(=O)O",                # L-arginine
        "N[C@@H](Cc1c[nH]cn1)C(=O)O"                # L-histidine
    ]
    
    # Precompute a set of canonical SMILES for standard proteinogenic amino acids
    prot_smiles_set = set()
    for aas in proteinogenic_aas:
        aas_mol = Chem.MolFromSmiles(aas)
        if aas_mol:
            canon = Chem.MolToSmiles(aas_mol, isomericSmiles=True)
            prot_smiles_set.add(canon)
    
    # For debugging one might want to see the input canonical SMILES.
    # print("Input canonical SMILES:", input_canonical)
    
    # If the input molecule exactly matches one of the standard amino acids, then it is proteinogenic.
    if input_canonical in prot_smiles_set:
        return False, "Matches a standard proteinogenic amino acid"
    
    # Otherwise, if an amino acid backbone is found and the structure is not the same as any
    # canonical proteinogenic amino acid, we classify it as non-proteinogenic.
    return True, "Has amino acid backbone and does not match any standard proteinogenic amino acid"

# For testing purposes (comment out if using as a module):
if __name__ == "__main__":
    test_smiles_list = [
        "N[C@@H](CC1=CC=C(F)C=C1)C(O)=O",   # 4-fluorophenyl-L-alanine
        "N[C@@H](CCCCN)C(O)=O"              # L-ornithine (should be proteinogenic)
    ]
    for smi in test_smiles_list:
        result, reason = is_non_proteinogenic_amino_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")