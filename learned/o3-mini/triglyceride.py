"""
Classifies: CHEBI:17855 triglyceride
"""
"""
Classifies: Triglycerides – any glyceride resulting from the condensation of all three hydroxy groups of glycerol (propane-1,2,3-triol) with fatty acids.
The function is_triglyceride(smiles: str) returns a boolean and a string explaining the result.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is defined as a glyceride formed by the condensation of all three hydroxy groups 
    of glycerol with fatty acids. In essence:
      - The molecule should be valid.
      - It should contain a glycerol backbone. Here we approximate the glycerol backbone as a 3‐carbon chain where
        the terminal carbons (CH2) and the middle carbon (CH) are present.
      - It must have exactly three ester bonds (-O-C(=O)-), which come from the esterification of glycerol’s HO groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triglyceride, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES to molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Heuristic check for a glycerol backbone.
    # We look for a chain of three carbons with the appropriate hybridization.
    # Note: Even after esterification, the CH2 and CH atoms remain.
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone (propane-1,2,3-triol core) not found"

    # Find ester groups. In a triglyceride each hydroxyl is esterified forming an ester bond -O-C(=O)-.
    # We use a SMARTS pattern to detect ester moieties.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Expected exactly 3 ester bonds, found {len(ester_matches)}"

    # For further confidence, we may check that the molecule has a high enough molecular weight 
    # (most triglycerides are > 500 Da) and a reasonable count of carbon atoms (at least 20).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical triglyceride"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Too few carbon atoms ({c_count}) for a typical triglyceride"

    return True, "Contains a glycerol backbone fully esterified with three fatty acid chains"

# Example usage (uncomment to test):
# smiles_list = [
#     "O(C(=O)CCCCCCCCCCCCC(C)C)[C@H](COC(=O)CCCCCCCCCCCCCCCCCCC)COC(=O)CCCCCCC",  # TG(8:0/i-16:0/20:0)
#     "CCCCCCCC\\C=C/CCCCCCCC(=O)OCC(COC(=O)CCCCCCC\\C=C/CCCCCCCC)OC(=O)CCCCCCC\\C=C/CCCCCCCC"  # triolein
# ]
# for smiles in smiles_list:
#     result, reason = is_triglyceride(smiles)
#     print(f"SMILES: {smiles}\nTriglyceride? {result} ({reason})\n")