"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
Definition: A heterocyclic fatty acid containing an epoxide ring as part of its structure.
A valid epoxy fatty acid must:
  - Have a free carboxylic acid group.
  - Contain an overall sufficient number of carbon atoms (we require at least 15).
  - Possess a long “linear” aliphatic chain (at least 12 carbons beyond the acid carbon)
    in which an epoxide ring (a three‐membered heterocycle formed by two carbons and an oxygen)
    is found at an internal (non‐terminal) position.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as an epoxy fatty acid based on its SMILES string.

    Strategy:
      1. Parse the molecule.
      2. Verify there is at least one free carboxylic acid group.
      3. Check that the overall molecule contains a minimum total carbon count.
      4. Identify epoxide rings (three-membered cycles with two carbons and one oxygen);
         record the carbon atoms involved.
      5. For each carboxylic acid group, start from its carbonyl carbon and for each attached 
         carbon that is "allowed" (i.e. an sp3-like carbon that is not embedded in a ring except if 
         it belongs to a 3-membered ring), perform a depth-first search (DFS) to enumerate all simple 
         chain paths.
      6. For each path (chain) found with at least 12 carbon atoms (not counting the acid carbon):
          - Check whether at least one internal (i.e. non-terminal) carbon on that chain is
            also part of an epoxide ring.
      7. If such a chain exists, the molecule qualifies as an epoxy fatty acid.
      
    Args:
        smiles (str): SMILES string representing the molecule.

    Returns:
        bool: True if the molecule qualifies as an epoxy fatty acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # (1) Parse the molecule from its SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (2) Check for at least one free carboxylic acid group.
    # We look for the acid group pattern: carbonyl attached to an –OH.
    acid_smarts = Chem.MolFromSmarts("C(=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    if not acid_matches:
        return False, "Missing carboxylic acid group (required fatty acid functionality)"

    # (3) Check overall number of carbons (must be at least 15).
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) < 15:
        return False, f"Too few carbon atoms ({len(carbons)}); must be at least 15 for a fatty acid"

    # (4) Identify epoxide rings.
    # Epoxide pattern: three connected atoms in a ring: carbon-oxygen-carbon, all in a 3-membered ring.
    epoxide_pattern = Chem.MolFromSmarts("[C;r3][O;r3][C;r3]")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if not epoxide_matches:
        return False, "No epoxide ring (three-membered heterocycle with oxygen) detected."
    
    # Record the indices of carbon atoms that are part of any epoxide ring.
    epoxide_carbons = set()
    for match in epoxide_matches:
        # match: (C, O, C) where the oxygen is at position 1.
        epoxide_carbons.add(match[0])
        epoxide_carbons.add(match[2])
    
    # Pre-calculate ring information for allowed atom checks.
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]

    def is_allowed_chain_atom(atom):
        """
        Returns True if an atom is a carbon that can be part of the linear aliphatic chain.
        Allowed if:
          - It is carbon (atomic number 6) AND
          - It is not in any ring larger than 3 atoms,
            OR if it is in a ring then every ring it falls in has exactly 3 atoms.
        """
        if atom.GetAtomicNum() != 6:
            return False
        idx = atom.GetIdx()
        if not atom.IsInRing():
            return True
        # Allow atom if it is only present in rings of size 3.
        for r in rings:
            if idx in r and len(r) > 3:
                return False
        return True

    # (5) Define DFS to find chain paths starting from a given atom.
    # We require that from the acid group neighbor, we can follow a contiguous chain (without repeating atoms)
    # that has at least chain_threshold carbons and contains at least one epoxide carbon in an internal position.
    chain_threshold = 12  # required number of chain carbons (excluding the acid carbon)

    def dfs_chain(current_idx, visited, path):
        """
        Depth-first search from the current atom.
        'path' is the list (chain) of atom indices traversed so far.
        'visited' is the set of visited atom indices.
        Returns True as soon as a qualifying chain is found.
        """
        # If the chain is sufficiently long, check whether an internal carbon is part of an epoxide.
        if len(path) >= chain_threshold:
            # Only consider internal atoms, i.e. excluding the first and last positions.
            internal_atoms = path[1:-1]
            if any(atom_idx in epoxide_carbons for atom_idx in internal_atoms):
                return True
        current_atom = mol.GetAtomWithIdx(current_idx)
        for nbr in current_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            if not is_allowed_chain_atom(nbr):
                continue
            visited.add(nbr_idx)
            if dfs_chain(nbr_idx, visited, path + [nbr_idx]):
                return True
            visited.remove(nbr_idx)
        return False

    # (6) For each acid match, try to find a valid chain.
    # The acid match provides a carboxylic acid group; we use the carbonyl carbon as the starting point.
    valid_chain_found = False
    for match in acid_matches:
        acid_carbon_idx = match[0]  # the carbonyl carbon (C in C(=O)[O;H])
        acid_atom = mol.GetAtomWithIdx(acid_carbon_idx)
        # Explore each neighbor of the acid carbon that is a carbon and allowed.
        for nbr in acid_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if not is_allowed_chain_atom(nbr):
                continue
            visited = {acid_carbon_idx, nbr.GetIdx()}
            # Start the DFS chain with this neighbor.
            if dfs_chain(nbr.GetIdx(), visited, [nbr.GetIdx()]):
                valid_chain_found = True
                break
        if valid_chain_found:
            break

    if not valid_chain_found:
        return False, ("Aliphatic chain too short or epoxide ring not embedded in a long, linear chain "
                       "from the acid group")
    
    return True, ("Molecule has a free carboxylic acid group and a long "
                  "aliphatic chain (≥12 carbons) containing an internal epoxide ring.")

# For testing: you can run this file directly to try one example.
if __name__ == "__main__":
    test_smiles = "O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H](CCCCC)O1)O)O"
    result, reason = is_epoxy_fatty_acid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)