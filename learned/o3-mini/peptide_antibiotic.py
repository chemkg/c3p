"""
Classifies: CHEBI:25903 peptide antibiotic
"""
#!/usr/bin/env python
"""
Classifies: Peptide Antibiotics
A chemically diverse class of peptides that exhibit antimicrobial properties.

This function uses refined heuristic rules to capture peptide‐like backbones while
minimizing false positives from small or simple peptides.
Refined criteria:
  1. The molecule must be successfully parsed.
  2. It must contain at least 3 amide bonds (C(=O)N) as a proxy for peptide bonds.
  3. It must have at least 2 peptide-backbone motifs 
     (chiral centers in the pattern "[C@H](N)C(=O)" or "[C@@H](N)C(=O)").
  4. Its molecular weight must be at least 500 Da (to filter out very short peptides).
  5. It must contain at least 4 nitrogen atoms (consistent with an amino acid–rich structure).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string using heuristics.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the peptide antibiotic heuristics, False otherwise.
        str: A message explaining the reason for the classification.
    """
    # 1. Attempt to parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Check for amide bonds (a proxy for peptide bonds).
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    count_amide = len(amide_matches)
    if count_amide < 3:
        return False, f"Insufficient amide bonds ({count_amide}); likely too short to be a peptide antibiotic"
    
    # 3. Look for peptide-backbone motifs.
    # These patterns detect a chiral alpha carbon attached to an amino group and a carbonyl.
    backbone_pattern1 = Chem.MolFromSmarts("[C@H](N)C(=O)")
    backbone_pattern2 = Chem.MolFromSmarts("[C@@H](N)C(=O)")
    matches1 = mol.GetSubstructMatches(backbone_pattern1)
    matches2 = mol.GetSubstructMatches(backbone_pattern2)
    count_backbone = len(matches1) + len(matches2)
    if count_backbone < 2:
        return False, f"Insufficient peptide backbone motifs ({count_backbone}); need at least 2 to suggest a multi-residue peptide"
    
    # 4. Check the molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical peptide antibiotic"
    
    # 5. Count the number of nitrogen atoms; peptides possess several amino groups.
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 4:
        return False, f"Too few nitrogen atoms ({nitrogen_count}); unexpected for a peptide backbone"
    
    return True, ("Contains multiple amide bonds and peptide backbone motifs, adequate molecular weight, "
                  "and sufficient nitrogen content consistent with peptide antibiotics")

# Example usage: (if running this module directly)
if __name__ == "__main__":
    # Test with a known peptide antibiotic SMILES (e.g., surfactin D)
    test_smiles = "[H][C@@]1(CCCCCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)O1"
    result, message = is_peptide_antibiotic(test_smiles)
    print(result, message)