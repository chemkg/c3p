"""
Classifies: CHEBI:16158 steroid sulfate
"""
#!/usr/bin/env python
"""
Classifies: steroid sulfate
Definition: A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid with sulfuric acid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    The molecule must contain:
      - A steroid nucleus (approximated by a four-ring fused system SMARTS).
      - A sulfate ester group (OS(=O)(=O)) attached to a carbon (i.e. via a former hydroxy group) on the steroid nucleus.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a steroid sulfate, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Define a SMARTS pattern for a steroid nucleus.
    # This is an approximate core for the tetracyclic steroid framework.
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C1CCC2C3")
    if steroid_pattern is None:
        return False, "Error in steroid SMARTS pattern."
        
    # Check if the molecule contains a steroid-like fused ring system.
    steroid_matches = mol.GetSubstructMatches(steroid_pattern)
    if not steroid_matches:
        return False, "No steroid-like fused ring system found."
        
    # For further checking, we take one steroid match (set of atom indexes)
    steroid_core_atoms = set(steroid_matches[0])
        
    # Define a SMARTS pattern for a sulfate ester group.
    # This looks for an oxygen bound to an S(=O)(=O) fragment.
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)")
    if sulfate_pattern is None:
        return False, "Error in sulfate SMARTS pattern."
        
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate ester group (OS(=O)(=O)) found."
        
    # Check that at least one sulfate group is attached via its oxygen to 
    # a carbon atom of the steroid core. This indicates that the sulfate is an ester of a steroid hydroxy.
    for match in sulfate_matches:
        # By our SMARTS, match[0] should be the oxygen that is directly attached to the S.
        sulfate_oxygen = mol.GetAtomWithIdx(match[0])
        # Look at its neighbors (besides the sulfate S)
        for nbr in sulfate_oxygen.GetNeighbors():
            # Exclude the sulfur (atomic number 16)
            if nbr.GetAtomicNum() == 16:
                continue
            # Check if this neighboring atom is a carbon and part of the steroid nucleus.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in steroid_core_atoms:
                return True, "Molecule contains a steroid nucleus with a sulfate ester attached to a hydroxy group."
                
    return False, "Sulfate ester not attached to steroid nucleus as required."

# For testing purposes: (you can run these lines to perform a simple test)
if __name__ == "__main__":
    test_smiles = "O(S(O)(=O)=O)[C@@H]1CC=2[C@]([C@]3(CC[C@]4([C@]([C@@]3(CC2)[H])(CC[C@@]4([C@H](C)CCCC(C)C)[H])[H])C)[H])(C)CC1"
    result, reason = is_steroid_sulfate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)