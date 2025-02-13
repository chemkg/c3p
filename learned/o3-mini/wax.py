"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: Wax
Definition: A wax is an organic compound or mixture of compounds composed of long‐chain molecules that are malleable at ambient temperatures.
This classifier applies several heuristics:
  - The molecule must yield a valid RDKit object.
  - It must only contain allowed elements: C, H, O.
  - Molecular weight must be at least 300 Da.
  - There must be at least 15 carbon atoms.
  - The fraction of aromatic carbons must be below 30%.
  - The fraction of carbons in rings must be below 50%.
  - There must be at least one uninterrupted carbon chain of at least 12 carbon atoms.
  - There must be at least 4 rotatable bonds.
  - Finally, to help distinguish wax esters from free acids or multiester lipids (like triacylglycerols), the molecule must contain exactly one ester bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from functools import lru_cache

def longest_carbon_chain_length(mol):
    """
    Computes the longest continuous chain (path) of only carbon atoms.
    Uses a DFS on a subgraph composed solely of carbon atoms, with memoization.

    Args:
        mol (rdkit.Chem.Mol): Molecule.

    Returns:
        int: The number of carbons in the longest uninterrupted carbon chain.
    """
    # Get indices of all carbon atoms.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_indices:
        return 0

    # Build an adjacency list for carbon atoms only.
    carbon_neighbors = {idx: [] for idx in carbon_indices}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors[idx].append(neighbor.GetIdx())

    @lru_cache(maxsize=None)
    def dfs(current, visited):
        max_length = 1
        for neighbor in carbon_neighbors[current]:
            if neighbor not in visited:
                new_length = 1 + dfs(neighbor, visited + (neighbor,))
                if new_length > max_length:
                    max_length = new_length
        return max_length

    longest = 0
    for start in carbon_indices:
        chain_length = dfs(start, (start,))
        if chain_length > longest:
            longest = chain_length
    return longest

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string by applying a set of heuristic rules:
      1. Valid molecule must be formed from allowed elements (C, H, O).
      2. The molecular weight must be at least 300 Da.
      3. It must contain at least 15 carbon atoms.
      4. The aromatic carbon fraction must be below 30%.
      5. The fraction of carbons in rings must be below 50%.
      6. There must be at least one uninterrupted carbon chain of 12 or more carbons.
      7. To avoid misclassifying fatty acids or multiester lipids, the molecule must contain exactly one ester bond,
         defined as the substructure [#6]-C(=O)-O-[#6].
      8. At least 4 rotatable bonds must be present (to indicate flexibility/malleability).

    Args:
        smiles (str): SMILES string representing the molecule.

    Returns:
        bool: True when the molecule meets all the wax criteria, False otherwise.
        str: A textual explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Allowed elements: Carbon, Hydrogen, Oxygen.
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element: {atom.GetSymbol()}."

    # Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low for a wax ({mol_wt:.1f} Da)."

    # Count carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbon_atoms)
    if num_carbons < 15:
        return False, f"Insufficient carbon count ({num_carbons}). Must have at least 15 carbons for a wax."

    # Calculate aromatic carbon fraction.
    aromatic_carbons = [atom for atom in carbon_atoms if atom.GetIsAromatic()]
    frac_aromatic = len(aromatic_carbons) / num_carbons
    if frac_aromatic >= 0.3:
        return False, f"Too many aromatic carbons ({len(aromatic_carbons)} out of {num_carbons})."

    # Calculate the fraction of carbons in rings.
    ring_carbons = [atom for atom in carbon_atoms if atom.IsInRing()]
    frac_ring = len(ring_carbons) / num_carbons
    if frac_ring >= 0.5:
        return False, f"Too many carbons are in rings ({len(ring_carbons)} out of {num_carbons})."

    # Check for long uninterrupted carbon chain.
    longest_chain = longest_carbon_chain_length(mol)
    if longest_chain < 12:
        return False, f"Longest continuous carbon chain is only {longest_chain} atoms (requires at least 12)."

    # Check number of rotatable bonds.
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds < 4:
        return False, f"Too few rotatable bonds ({rot_bonds}); molecule may be too rigid."

    # Count ester bonds.
    # We define an ester bond as a fragment with the pattern: [#6]-C(=O)-O-[#6]
    ester_pattern = Chem.MolFromSmarts("[#6]-C(=O)-O-[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected exactly one ester bond for a wax, found {len(ester_matches)}."

    return True, (f"Valid wax: Molecular weight {mol_wt:.1f} Da, {num_carbons} carbons "
                  f"(aromatic fraction {frac_aromatic:.2f}, ring fraction {frac_ring:.2f}), "
                  f"longest chain of {longest_chain} carbons, {rot_bonds} rotatable bonds, "
                  "and one ester bond detected.")

# Example usage (uncomment to test):
# examples = [
#     "O(CCCCCCCC/C=C\\C/C=C\\CCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\CCCCC",  # Linoleyl linoleate (TP)
#     "CCCCCCCCCCCCCCCC(=O)OCCCCCCCCC",  # nonyl palmitate (TP)
#     "OC(=O)CC\\C=C\\C\\C=C\\C\\C=C\\C\\C=C\\CCCCCCCC",  # Docosatetraenoic acid (FP)
# ]
# for s in examples:
#     valid, reason = is_wax(s)
#     print(s)
#     print(valid, reason)