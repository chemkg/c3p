"""
Classifies: CHEBI:25903 peptide antibiotic
"""
#!/usr/bin/env python
"""
Classifies: Peptide Antibiotics
A chemically diverse class of peptides that exhibit antimicrobial properties.
This function uses several heuristic rules to infer if a test molecule might belong to this class.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are a diverse set of peptides with antimicrobial properties.
    
    Heuristic criteria applied:
      1. The molecule must be successfully parsed.
      2. It must contain multiple amide bonds (C(=O)N), which are used as proxies for peptide bonds.
         (Peptides of biological relevance typically have at least 3 peptide bonds.)
      3. The molecular weight should be above a threshold (e.g., 300 Da) expected for peptides.
      4. The molecule should contain multiple nitrogen atoms, characteristic of amino acid backbones.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule meets the peptide antibiotic heuristics, False otherwise.
        str: A message explaining the reason for the result.
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Check for multiple amide bonds; peptide bonds show up as C(=O)N.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    count_amide = len(amide_matches)
    if count_amide < 3:
        return False, f"Insufficient amide bonds ({count_amide}) to suggest a peptide chain"

    # 3. Check molecular weight; many peptide antibiotics are >300 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical peptide antibiotic"

    # 4. Check the number of nitrogen atoms (many amino acids add nitrogen).
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 4:
        return False, f"Too few nitrogen atoms ({nitrogen_count}) to represent a peptide backbone"
    
    # All heuristics passed – this molecule is likely to be a peptide antibiotic.
    return True, "Contains multiple amide bonds, adequate molecular weight, and sufficient nitrogen atoms consistent with peptide antibiotics"