"""
Classifies: CHEBI:18000 aralkylamine
"""
#!/usr/bin/env python3
"""
Classifies: aralkylamine
Definition: An alkylamine in which the alkyl group is substituted by an aromatic group.
This program checks if a molecule (given as a SMILES string) contains an aliphatic amine N
that is substituted with a carbon chain which eventually leads to an aromatic (ring) carbon.
"""

from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine having an alkyl substituent that bears an aromatic group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True if the molecule is an aralkylamine, False otherwise,
                     with a reason for the classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper: limited-depth DFS that only traverses aliphatic carbons (non-aromatic)
    # and looks for an aromatic carbon.
    def dfs_find_aromatic(atom, depth, visited):
        if depth < 0:
            return False
        visited.add(atom.GetIdx())
        # Check all neighbors
        for nbr in atom.GetNeighbors():
            # We do not want to return immediately back to the amine nitrogen; but in our DFS, visited will help.
            if nbr.GetIdx() in visited:
                continue
            # If we find a carbon that is aromatic, we have reached an aromatic ring.
            if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic():
                return True
            # Continue along non-aromatic carbons.
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                if dfs_find_aromatic(nbr, depth - 1, visited):
                    return True
        return False

    # Look for nitrogens in the molecule that are not aromatic (i.e. part of an alkylamine)
    for atom in mol.GetAtoms():
        # Consider nitrogen atoms that are not aromatic.
        if atom.GetAtomicNum() == 7 and not atom.GetIsAromatic():
            # For an aralkylamine the substituent attached to N must be an aliphatic carbon (i.e. not aromatic).
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                    # We allow a chain length up to 3 bonds. Copy the DFS (with new visited set) starting from this neighbor.
                    if dfs_find_aromatic(nbr, depth=3, visited=set()):
                        return True, "Contains an aliphatic amine nitrogen with a substituent that connects to an aromatic ring"
    
    return False, "No aralkyl substituent detected on an alkylamine nitrogen"

# For testing purposes you can run:
if __name__ == "__main__":
    # Example: benzylamine is the simplest aralkylamine (SMILES: "NCc1ccccc1")
    test_smiles = "NCc1ccccc1"
    result, reason = is_aralkylamine(test_smiles)
    print("SMILES:", test_smiles)
    print("Aralkylamine?", result)
    print("Reason:", reason)