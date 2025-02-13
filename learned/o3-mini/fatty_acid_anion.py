"""
Classifies: CHEBI:28868 fatty acid anion
"""
#!/usr/bin/env python3
"""
Classifies: Fatty acid anion
Definition: The conjugate base of a fatty acid, arising from deprotonation 
of the carboxylic acid group of the corresponding fatty acid.
A fatty acid anion is expected to have one terminal carboxylate group ([C(=O)[O-]])
attached to an aliphatic carbon chain of reasonable length.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    
    Heuristic criteria:
      - The molecule must be parsed correctly.
      - It must contain exactly one deprotonated carboxylate group ([CX3](=O)[O-]).
      - The carboxylate group must be terminal, that is, its carbon is attached 
        to exactly one carbon neighbor.
      - The carbon chain attached (alpha carbon) must have a reasonable length.
      - The molecule overall should have a modest number of C atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a fatty acid anion, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a carboxylate group.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # The fatty acid anion should have exactly one carboxylate group.
    if len(carboxylate_matches) == 0:
        return False, "No deprotonated carboxylate group ([C(=O)[O-]]) found"
    if len(carboxylate_matches) > 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups; expected exactly one"

    # Get the matched indices: by our SMARTS definition the first atom is the carbon of the carboxylate.
    match = carboxylate_matches[0]
    carboxyl_carbon_idx = match[0]
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

    # Check that the carboxylate carbon is terminal (should have exactly one carbon neighbor).
    carbon_neighbors = [nbr for nbr in carboxyl_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylate group is not terminal (expected exactly one carbon neighbor)"

    # For a minimal fatty acid, require that there are enough carbon atoms.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 5:
        return False, f"Too few carbon atoms ({total_carbons}) to be a fatty acid"

    # From the alpha-carbon (the carbon neighbor attached to the carboxylate group)
    # we try to measure the chain length using a simple DFS that only traverses carbon atoms.
    alpha_carbon = carbon_neighbors[0]

    # Helper: recursively compute the maximum chain length (number of carbon atoms) from a starting carbon.
    def get_longest_chain(atom, visited):
        visited.add(atom.GetIdx())
        max_length = 1  # count current atom
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                # Avoid stepping back into the carboxylate carbon.
                chain_length = 1 + get_longest_chain(nbr, visited.copy())
                if chain_length > max_length:
                    max_length = chain_length
        return max_length

    chain_length = get_longest_chain(alpha_carbon, set())
    # Set a minimal chain length threshold. Here we choose 4 as a minimal number of carbons 
    # in the chain (not including the carboxylate carbon) though many fatty acids are longer.
    if chain_length < 4:
        return False, f"Aliphatic chain appears too short (chain length = {chain_length})"

    # Optionally, check that the molecular weight is within a range typical for fatty acids 
    # (most fatty acids have moderate molecular weights). Here we allow a wide range so we skip this.

    return True, "Contains a terminal carboxylate group attached to a sufficient aliphatic chain, consistent with a fatty acid anion"

# Example usage (if run as a script, one can test with one of the provided SMILES strings):
if __name__ == "__main__":
    # Example SMILES: hexadecanoate
    test_smiles = "CCCCCCCCCCCCCCCC([O-])=O"
    is_fa, reason = is_fatty_acid_anion(test_smiles)
    print(is_fa, reason)