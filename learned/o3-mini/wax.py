"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: Wax
Definition: A wax is an organic compound or mixture of compounds composed of long‐chain molecules and is malleable at ambient temperatures.
This classifier uses an improved set of heuristics:
  - The molecule must parse correctly.
  - The molecular weight must be above 300 Da.
  - The molecule must have at least 15 carbon atoms.
  - The molecule must be mostly aliphatic: the fraction of aromatic carbons must be below 30%
    and the fraction of carbons located in any rings must be below 50%.
  - The molecule must contain only allowed elements (C, H, O).
  - Contains at least one uninterrupted carbon chain (pure C–C bonds) of 12 or more carbons.
  - Must have at least 4 rotatable bonds (to help indicate “malleability”).
  
These updated heuristics help avoid classification of highly cyclic (and/or heteroatom‐rich) structures.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from functools import lru_cache

def longest_carbon_chain_length(mol):
    """
    Computes the longest continuous path containing only carbon atoms.
    A DFS is used on the subgraph of carbon atoms only. We use memoization for efficiency.
    
    Args:
        mol (rdkit.Chem.Mol): Molecule object.

    Returns:
        int: Length (number of carbons) of the longest uninterrupted carbon chain.
    """
    # Get the indices of all carbon atoms.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_indices:
        return 0

    # Build an adjacency list (only carbons) from the molecule.
    carbon_neighbors = {idx: [] for idx in carbon_indices}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors[idx].append(neighbor.GetIdx())

    @lru_cache(maxsize=None)
    def dfs(current, visited):
        max_length = 1  # count current atom
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
    Determines if a molecule is a wax based on its SMILES string by applying a suite
    of heuristic rules:
      1. The input must yield a valid molecule.
      2. Molecular weight must be at least 300 Da.
      3. The molecule must contain at least 15 carbon atoms.
      4. The aromatic carbon fraction must be below 30%.
      5. The fraction of carbons in rings must be below 50%.
      6. The molecule should contain only allowed elements: C, H, and O.
      7. There must be an uninterrupted carbon chain of at least 12 atoms.
      8. The molecule should be flexible, with at least 4 rotatable bonds.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a wax, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check for forbidden elements: allow only C (6), O (8) and implicit H.
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element: {atom.GetSymbol()}."

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low for a wax ({mol_wt:.1f} Da)."
    
    # Count carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbon_atoms)
    if num_carbons < 15:
        return False, f"Insufficient carbon count ({num_carbons}). Must have at least 15 carbons for a wax."
    
    # Calculate aromatic carbon fraction
    aromatic_carbons = [atom for atom in carbon_atoms if atom.GetIsAromatic()]
    frac_aromatic = len(aromatic_carbons) / num_carbons
    if frac_aromatic >= 0.3:
        return False, f"Too many aromatic carbons ({len(aromatic_carbons)} out of {num_carbons})."
    
    # Calculate fraction of carbon atoms in rings
    ring_carbons = [atom for atom in carbon_atoms if atom.IsInRing()]
    frac_ring = len(ring_carbons) / num_carbons
    if frac_ring >= 0.5:
        return False, f"Too many carbons are in rings ({len(ring_carbons)} out of {num_carbons})."
    
    # Check for a long continuous carbon chain (must be at least 12)
    longest_chain = longest_carbon_chain_length(mol)
    if longest_chain < 12:
        return False, f"Longest continuous carbon chain is only {longest_chain} atoms (requires at least 12)."
    
    # Check number of rotatable bonds (to imply flexibility/malleability)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds < 4:
        return False, f"Too few rotatable bonds ({rot_bonds}); molecule may be too rigid."
    
    return True, f"Valid wax: Molecular weight {mol_wt:.1f} Da, {num_carbons} carbons (aromatic fraction {frac_aromatic:.2f}, ring fraction {frac_ring:.2f}), longest chain of {longest_chain} carbons, {rot_bonds} rotatable bonds."

# Example usage (uncomment the lines below to test):
# examples = [
#     "O(CCCCCCCC/C=C\\C/C=C\\CCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\CCCCC",  # Linoleyl linoleate (TP)
#     "CCCCCCCCCCCCCCCC(=O)OCCCCCCCCC",  # nonyl palmitate (TP)
#     "O=C(CCCCCCCCC)C=1C=CC(=NC1)CCCCCCCCC",  # 1-(6-Nonylpyridin-3-yl)decan-1-one (FP: contains N)
# ]
# for s in examples:
#     result, reason = is_wax(s)
#     print(s)
#     print(result, reason)