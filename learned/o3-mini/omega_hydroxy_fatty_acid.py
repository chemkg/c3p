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

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    The criteria:
      - It has a carboxylic acid group (-C(=O)O).
      - It has at least one terminal (end-of-chain) –OH group that is not part of the carboxyl group.
    
    The approach is:
      1. Parse the SMILES and add explicit hydrogens.
      2. Find carboxylic acid groups using SMARTS (which detects -C(=O)[O;H]).
         Record the acid carbon(s).
      3. Find carbon atoms that are terminal (only one connected carbon neighbor).
      4. Exclude the terminal carbon that is part of the acid group.
      5. Check if the remaining terminal carbon(s) have a directly attached hydroxyl group (-OH)
         (a single-bonded oxygen with at least one hydrogen) that is not already assigned to the acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if molecule meets the criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to the molecule so that -OH hydrogens are represented.
    mol = Chem.AddHs(mol)

    # Identify the carboxylic acid group.
    # The SMARTS “C(=O)[O;H]” will match the carboxyl carbon and its hydroxyl oxygen.
    acid_smarts = "C(=O)[O;H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Missing carboxylic acid group (-C(=O)O)"
    
    # Record all acid carbon indices (the first atom of the match is the carboxyl carbon).
    acid_carbon_idxs = set(match[0] for match in acid_matches)
    # Also record the oxygen atoms involved in the acid group.
    acid_oxygen_idxs = set()
    for match in acid_matches:
        # The pattern is expected to yield (carbon, oxygen) or (carbon, carbonyl_oxygen, hydroxyl oxygen) depending on the pattern.
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1:
                acid_oxygen_idxs.add(idx)

    # Now, identify all carbon atoms (atomic number 6) that are terminal.
    # A terminal carbon (in a linear fatty acid chain) should be connected to exactly one other carbon.
    terminal_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Count the number of neighboring carbon atoms.
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_carbons.append(atom.GetIdx())
    
    # Remove the acid carbon from the terminal candidates.
    candidate_terminal_idxs = [idx for idx in terminal_carbons if idx not in acid_carbon_idxs]
    
    if not candidate_terminal_idxs:
        return False, "No terminal carbon (other than the carboxyl carbon) found in the chain"

    # For each candidate terminal carbon, check if it has an attached hydroxyl (-OH) that is not part of the acid group.
    omega_hydroxyl_found = False
    for c_idx in candidate_terminal_idxs:
        carbon = mol.GetAtomWithIdx(c_idx)
        for nbr in carbon.GetNeighbors():
            # We look for oxygen neighbors.
            if nbr.GetAtomicNum() != 8:
                continue
            # Skip if this oxygen is part of the acid group.
            if nbr.GetIdx() in acid_oxygen_idxs:
                continue
            # Check the bond type is single.
            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Make sure oxygen looks like an alcohol (has at least one hydrogen).
            if nbr.GetTotalNumHs() < 1:
                continue
            omega_hydroxyl_found = True
            break
        if omega_hydroxyl_found:
            break

    if not omega_hydroxyl_found:
        return False, "No terminal hydroxyl (-OH) group found on the opposite end of the chain from the carboxyl group"

    return True, ("Found a carboxylic acid (-C(=O)O) group and a terminal hydroxyl (-OH) group on a terminal carbon "
                  "distinct from the acid carbon, consistent with an omega-hydroxy fatty acid.")

# Example usage (uncomment for testing):
# if __name__ == "__main__":
#     test_smiles = "OCCCCCCCC(O)=O"  # Example: 8-hydroxyoctanoic acid
#     result, reason = is_omega_hydroxy_fatty_acid(test_smiles)
#     print(result, reason)