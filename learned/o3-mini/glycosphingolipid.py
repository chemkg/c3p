"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: Glycosphingolipid
A glycosphingolipid is a glycolipid that is a carbohydrate-containing derivative of a sphingoid or ceramide.
The carbohydrate residue is attached by a glycosidic linkage to O-1 of the sphingoid.
This routine uses heuristic SMARTS patterns to look for a pyranose ring (sugar), an amide group typical of ceramides,
and evidence for a glycosidic linkage (via a ring oxygen connecting to a non-ring atom).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    Heuristics used:
      - The molecule must be a valid structure.
      - It must contain at least one pyranose ring (a six-membered ring with one oxygen).
      - It must contain an amide group (NC(=O)) typical of ceramides.
      - At least one oxygen in the sugar should be linked to an atom outside the ring (a crude check for a glycosidic linkage).
      - There should be a sizable number of aliphatic carbons typical of a lipid tail.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if it satisfies our criteria for a glycosphingolipid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Heuristic 1: Look for at least one pyranose sugar ring.
    # We use a flexible SMARTS for a saturated 6-membered ring that contains exactly one oxygen.
    # [OX2r6] ensures an sp3 oxygen in a ring of size 6 and [CX4r6] matches sp3 carbons in the ring.
    sugar_smarts = "[OX2r6]1[CX4r6][CX4r6][CX4r6][CX4r6][CX4r6]1"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No pyranose sugar ring found"
        
    # Heuristic 2: Look for an amide group typical of ceramides.
    # We use a simple SMARTS to match "N-C(=O)".
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("N[C;!$(C(=O))]C(=O)") ) and not mol.HasSubstructMatch(Chem.MolFromSmarts("NC(=O)")):
        # Using the relaxed pattern "NC(=O)" first.
        return False, "No amide group (typical of ceramide) found"
    # Alternatively you can check with a stricter pattern:
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("NC(=O)")):
        return False, "Amide group not found"

    # Heuristic 3: Check for a glycosidic linkage.
    # In a glycosphingolipid, one oxygen of the sugar ring should be linked (via a single bond)
    # to an atom that is not part of the sugar ring.
    glyco_linkage_found = False
    for match in sugar_matches:
        sugar_atom_indices = set(match)
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            # Check for oxygen atoms in the ring
            if atom.GetSymbol() == "O":
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in sugar_atom_indices:
                        glyco_linkage_found = True
                        break
            if glyco_linkage_found:
                break
        if glyco_linkage_found:
            break
    if not glyco_linkage_found:
        return False, "No oxygen outside the sugar ring detected as a glycosidic linkage"
        
    # Optional: Check that there is a sufficiently long aliphatic chain (lipid tail)
    aliphatic_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
    if len(aliphatic_carbons) < 10:
        return False, "Not enough aliphatic carbons to indicate a lipid tail"
        
    # If all heuristics pass, classify as a glycosphingolipid.
    return True, "Contains a pyranose sugar ring with a glycosidic linkage and an amide group (ceramide-like) plus a lipid tail"

# Example usage (can be removed or commented out when used as a module):
if __name__ == "__main__":
    # Try one of the example SMILES strings:
    smiles_example = "CCCCCCCCCCCCCC[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)C(O)CCCC"
    result, reason = is_glycosphingolipid(smiles_example)
    print("Result:", result)
    print("Reason:", reason)