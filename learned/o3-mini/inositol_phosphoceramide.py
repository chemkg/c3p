"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: inositol phosphoceramide
Definition:
  A phosphosphingolipid in which an inositol residue and the ceramide moiety are linked via a phosphodiester bridge.
  The ceramide moiety contains substituents R(1) and R(2) which vary with different sphingoid bases and fatty acyl moieties.
  
We approximate detection by:
 - Looking for an inositol ring substructure,
 - Looking for a phosphate (phosphodiester) group, and
 - Identifying a ceramide substructure represented by an amide group (C(=O)N) attached to at least one long aliphatic chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    Criteria used (approximate):
      1. Presence of an inositol ring (six-membered ring with multiple hydroxyl groups)
      2. Phosphodiester bridge (a phosphate connected via ester bonds; we search for OP(=O)(O)O)
      3. Ceramide-like moiety: an amide group (C(=O)N) and a long aliphatic chain (assessed by total carbon count)
     
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for inositol ring.
    # This SMARTS approximates a myo-inositol ring pattern.
    inositol_smarts = "O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring detected"
    
    # Check for phosphodiester bridge.
    # We look for a phosphate group with three oxygens (one typically linking to the inositol).
    phosphate_smarts = "[O;D2]-P(=O)([O;D1])[O;D2]"  # this pattern looks for a phosphodiester moiety (O-P(=O)-O)
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphodiester bridge (phosphate group) detected"
    
    # Look for ceramide-like moiety: an amide group.
    amide_smarts = "C(=O)N"
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond detected to suggest a ceramide moiety"
    
    # Check for long aliphatic chain (fatty acyl chains) typical in ceramides.
    # We count the number of carbon atoms in the molecule.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 20:
        return False, "Total carbon count too low to suggest long fatty acyl chains in ceramide"
    
    # Optionally also check for rotatable bonds as a proxy for flexible, long chains.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Not enough rotatable bonds; fatty acyl chains may be too short"
    
    # Passed all checks
    return True, "Molecule contains an inositol ring, a phosphodiester bridge, and a ceramide-like amide bond with long aliphatic chains"

# Example usage (if run as main, uncomment):
# if __name__ == "__main__":
#     test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)C(O)CCCCCCCCCCCCCC"
#     result, reason = is_inositol_phosphoceramide(test_smiles)
#     print(result, reason)