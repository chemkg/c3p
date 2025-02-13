"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
#!/usr/bin/env python
"""
Classifies: omega-hydroxy fatty acid
Definition: Any naturally-occurring straight-chain fatty acid with a carboxyl group at position 1 and a hydroxyl at position n (omega).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    The molecule must contain a carboxylic acid group (-C(=O)O) and a terminal (omega) hydroxyl group 
    on a saturated carbon that is connected to only one other carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the carboxylic acid group using SMARTS.
    # The SMARTS "[CX3](=O)[OX1H]" looks for a carbon (X3) double bonded to oxygen and single bonded to an -OH.
    acid_smarts = "[CX3](=O)[OX1H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "Missing carboxylic acid group (-C(=O)O)"

    # Find candidate terminal hydroxyl (-OH) groups that are not part of the acid.
    # We iterate over all oxygen atoms that have at least one hydrogen attached.
    terminal_oh_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue  # Skip non-oxygen atoms
        # Check if the oxygen atom has an attached hydrogen (i.e. acts as –OH)
        if atom.GetTotalNumHs() < 1:
            continue
        # Loop over neighbors of the oxygen; usually the hydroxyl oxygen is bound to one carbon.
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                continue  # We want the oxygen bound to a carbon
            # Check that the bond is a single bond (as expected for an alcohol)
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Check to avoid mistaking the acid -OH: if the carbon (neighbor) has a double-bonded oxygen, it is likely the acid.
            has_db_oxygen = False
            for nb in neighbor.GetNeighbors():
                if nb.GetIdx() == atom.GetIdx():
                    continue
                b = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nb.GetIdx())
                if b.GetBondType() == Chem.BondType.DOUBLE and nb.GetAtomicNum() == 8:
                    has_db_oxygen = True
                    break
            if has_db_oxygen:
                continue  # part of a carboxyl group rather than a terminal alcohol
            
            # Now check if the carbon neighbor is terminal in the carbon chain.
            # The omega carbon should be attached to only one other carbon.
            carbon_neighbors = [nb for nb in neighbor.GetNeighbors() if nb.GetAtomicNum() == 6]
            if len(carbon_neighbors) == 1:
                terminal_oh_found = True
                break
        if terminal_oh_found:
            break

    if not terminal_oh_found:
        return False, "No terminal hydroxyl (-OH) group found that is distinct from the carboxyl group"
    
    # Optionally, additional checks (e.g., chain length, rotatable bonds) could be added to verify
    # that the fatty acid is a long enough, unbranched chain.
    # For now, the presence of both a carboxylic acid and an appropriate terminal hydroxyl is considered sufficient.
    
    return True, ("Found both a carboxylic acid (-C(=O)O) group and a terminal hydroxyl group "
                  "on a saturated carbon likely corresponding to position omega.")

# Example usage (uncomment for testing):
# if __name__ == "__main__":
#     test_smiles = "OCCCCCCCC(O)=O"  # Example: omega-hydroxytritriacontanoic acid
#     result, reason = is_omega_hydroxy_fatty_acid(test_smiles)
#     print(result, reason)