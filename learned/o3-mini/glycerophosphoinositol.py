"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: glycerophosphoinositol 
Definition: Any glycerophospholipid having the polar alcohol inositol 
esterified to the phosphate group at the sn-3 position of the glycerol backbone.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol has a glycerol backbone, a phosphate group and an inositol ring attached
    to the phosphate (i.e. forming a phosphoinositol headgroup).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Explanation for the classification decision
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a simplified SMARTS for the inositol ring (myo-inositol).
    # This pattern looks for a six-membered aliphatic ring with hydroxyl groups on each carbon.
    inositol_smarts = "C1C(O)C(O)C(O)C(O)C1O"
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if inositol_pattern is None:
        return False, "Error in the inositol SMARTS pattern"

    # Define a SMARTS for a phosphate group attached to an oxygen.
    # We do not require full stereochemistry here. We require that one of the substituents 
    # of phosphorus is an oxygen bound to an inositol ring.
    phosphate_smarts = "P(=O)(O)O"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if phosphate_pattern is None:
        return False, "Error in the phosphate SMARTS pattern"
    
    # Define a SMARTS representing a glycerol backbone.
    # This simplified pattern looks for two adjacent hydroxylated carbons ending in a CH2-O unit.
    glycerol_smarts = "C(O)C(O)CO"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    if glycerol_pattern is None:
        return False, "Error in the glycerol SMARTS pattern"

    # Check for the inositol ring presence.
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring found"

    # Check for a phosphate group.
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Check if at least one phosphate group is connected to an oxygen that is part of an inositol ring.
    # We loop over each phosphorus match and inspect its neighbor atoms.
    phosphoinositol_found = False
    for match in phosphate_matches:
        # match is a tuple with indices corresponding to atoms in the phosphate pattern.
        P_idx = match[0]  # phosphorus atom should be the first atom in our phosphate pattern SMARTS.
        P_atom = mol.GetAtomWithIdx(P_idx)
        # Look at neighbors of phosphorus; we expect one or more oxygen atoms.
        for neighbor in P_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # oxygen
                # Check if this oxygen (or the substructure containing it) is part of an inositol ring.
                # Create a one-atom query for this oxygen.
                # We use a temporary mol with just this atom to see if it is matched within an inositol substructure.
                # Alternatively, we check if this neighbor's index appears in any inositol match.
                for inositol_match in inositol_matches:
                    if neighbor.GetIdx() in inositol_match:
                        phosphoinositol_found = True
                        break
                if phosphoinositol_found:
                    break
        if phosphoinositol_found:
            break
    if not phosphoinositol_found:
        return False, "Phosphate group is not attached to an inositol ring"
    
    # Finally, check for a glycerol backbone (which is required for a glycerophospholipid).
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"
    
    return True, "Molecule contains a glycerol backbone with a phosphate group bearing an inositol headgroup"

# Example usage:
if __name__ == "__main__":
    test_smiles = "C(COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O)OC(=O)CCCCCCCCCCCCCCC"  # Simplified example
    result, reason = is_glycerophosphoinositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)