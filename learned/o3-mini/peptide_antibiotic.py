"""
Classifies: CHEBI:25903 peptide antibiotic
"""
#!/usr/bin/env python
"""
Classifies: Peptide Antibiotics
A chemically diverse class of peptides that exhibit antimicrobial properties.
This function uses several refined heuristic rules to better capture a peptide backbone.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    
    Refined heuristic criteria:
      1. The molecule must be successfully parsed.
      2. It must contain at least three amide bonds (C(=O)N) to indicate possible peptide bonds.
      3. It must show at least one peptide-backbone motif, for example "[C@H](N)C(=O)" or "[C@@H](N)C(=O)".
      4. Its molecular weight is expected to be above a set threshold (here, 400 Da).
      5. It must contain a sufficient number of nitrogen atoms (at least 4) as found in peptide backbones.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the peptide antibiotic heuristics, False otherwise.
        str: A message explaining the reason for the classification.
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Check for amide bonds (C(=O)N) as a general proxy for peptide bonds.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    count_amide = len(amide_matches)
    if count_amide < 3:
        return False, f"Insufficient amide bonds ({count_amide}) to suggest a peptide chain"
    
    # 3. Look for a peptide backbone pattern.
    #    These SMARTS look for a chiral alpha carbon bound to an amino group and connected to a carbonyl.
    backbone_pattern1 = Chem.MolFromSmarts("[C@H](N)C(=O)")
    backbone_pattern2 = Chem.MolFromSmarts("[C@@H](N)C(=O)")
    count_backbone = len(mol.GetSubstructMatches(backbone_pattern1)) + len(mol.GetSubstructMatches(backbone_pattern2))
    if count_backbone < 1:
        return False, "No recognizable peptide backbone pattern found"

    # 4. Check the molecular weight; many peptide antibiotics are relatively large.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical peptide antibiotic"
    
    # 5. Count the number of nitrogen atoms (an indicator of amino acid content).
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 4:
        return False, f"Too few nitrogen atoms ({nitrogen_count}) to represent a peptide backbone"
    
    return True, "Contains peptide backbone patterns, multiple amide bonds, adequate molecular weight, and sufficient nitrogen atoms consistent with peptide antibiotics"

# Example usage (if run as a script):
if __name__ == "__main__":
    # Use one of the true positive SMILES as a test (e.g., surfactin D).
    test_smiles = "[H][C@@]1(CCCCCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)O1"
    result, message = is_peptide_antibiotic(test_smiles)
    print(result, message)