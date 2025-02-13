"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
#!/usr/bin/env python
"""
Classifies: omega-hydroxy fatty acid
Definition: Any member of the class of naturally-occurring straight-chain fatty acids n carbon atoms long 
with a carboxyl group at position 1 and a hydroxyl at position n (omega).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    The molecule must contain a carboxylic acid group (-C(=O)O) and a terminal (omega) hydroxyl group 
    on a saturated carbon connected to only one other carbon.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens so that functional groups (e.g. -OH) are represented.
    mol = Chem.AddHs(mol)

    # Define a SMARTS pattern to detect the carboxylic acid group.
    # This pattern matches a carbonyl (C=O) attached to an -OH.
    acid_smarts = "C(=O)[O;H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Missing carboxylic acid group (-C(=O)O)"
    
    # Retrieve indices of oxygen atoms that are part of the carboxylic acid.
    acid_oh_idxs = set()
    for match in acid_matches:
        # In our SMARTS "C(=O)[O;H]", we expect the match to return a tuple where one of the atoms is the -OH oxygen.
        # We add any oxygen (atomic number 8) in the match that has at least one hydrogen.
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1:
                acid_oh_idxs.add(idx)

    # Look for candidate terminal hydroxyl (-OH) groups that are not part of the carboxylic acid.
    terminal_oh_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue  # Only consider oxygen atoms.
        # Skip if this oxygen is part of the carboxylic acid group.
        if atom.GetIdx() in acid_oh_idxs:
            continue
        # Ensure the oxygen has at least one hydrogen (i.e. acts as an -OH).
        if atom.GetTotalNumHs() < 1:
            continue
        # Examine neighbors, typically the hydroxyl oxygen is bound to one carbon.
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                continue  # We require the -OH to be attached to a carbon
            # Verify that the bond is single (as is expected for an alcohol)
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue

            # Check if this neighbor carbon is terminal in the chain.
            # Count the number of carbon atoms attached to this carbon.
            carbon_neighbors = [nb for nb in neighbor.GetNeighbors() if nb.GetAtomicNum() == 6]
            if len(carbon_neighbors) == 1:
                terminal_oh_found = True
                break
        if terminal_oh_found:
            break

    if not terminal_oh_found:
        return False, "No terminal hydroxyl (-OH) group found that is distinct from the carboxyl group"

    # Optionally, additional checks like chain length or unbranched structure could be added here.
    return True, ("Found both a carboxylic acid (-C(=O)O) group and a terminal hydroxyl (-OH) group "
                  "on a carbon with only one carbon neighbor, consistent with an omega-hydroxy fatty acid.")

# Example usage (uncomment for testing):
# if __name__ == "__main__":
#     test_smiles = "OCCCCCCCC(O)=O"  # Example: 8-hydroxyoctanoic acid
#     result, reason = is_omega_hydroxy_fatty_acid(test_smiles)
#     print(result, reason)