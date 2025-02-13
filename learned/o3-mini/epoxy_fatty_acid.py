"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
Definition: A heterocyclic fatty acid containing an epoxide ring as part of its structure.
A valid epoxy fatty acid must:
  - Have a free carboxylic acid group.
  - Contain an overall sufficient number of carbons (we require at least 15).
  - Have a long “linear” aliphatic chain (at least 12 carbons away from the acid carbon)
    that includes an epoxide (three‐membered heterocycle) in an internal (non‐terminal) position.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.

    The strategy is as follows:
      1. Parse the SMILES.
      2. Check that a free carboxylic acid group is present (using a SMARTS).
      3. Check that the total number of carbons is above a set threshold (>=15).
      4. Identify at least one epoxide ring (three-membered cycle consisting of two carbons and an oxygen)
         and record its carbon atoms.
      5. From the acid carbon (found from the acid SMARTS), perform a DFS along allowed chain carbons.
         Allowed chain carbons are those that are either not in any ring or are only in 3-membered rings.
      6. Collect all simple paths (no revisiting atoms) and determine if any path, representing the fatty acid
         chain (excluding the acid carbon), has at least 12 carbons AND includes at least one epoxide carbon,
         not at the very start or very end of the chain.
      7. If such a path exists then the molecule qualifies; otherwise it is rejected.

    Args:
        smiles (str): SMILES string representing the molecule.

    Returns:
        bool: True if the molecule qualifies as an epoxy fatty acid, False otherwise.
        str: Reason for classification.
    """
    # (1) Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (2) Check for the carboxylic acid group.
    acid_smarts = Chem.MolFromSmarts("C(=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    if not acid_matches:
        return False, "Missing carboxylic acid group (fatty acid functionality)"
    
    # (3) Check overall number of carbon atoms.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) < 15:
        return False, f"Too few carbon atoms ({len(carbons)}); not long enough to be a fatty acid"
    
    # (4) Identify epoxide ring(s): use the three-membered cycle of two carbons and one oxygen.
    epoxide_pattern = Chem.MolFromSmarts("[C;r3][O;r3][C;r3]")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if not epoxide_matches:
        return False, "No epoxide ring (three-membered heterocycle containing oxygen) detected."
    
    # Build a set of carbon atom indices that are part of any epoxide ring.
    epoxide_carbons = set()
    for match in epoxide_matches:
        # match is a tuple of three atoms [C, O, C] (the oxygen is match[1])
        epoxide_carbons.add(match[0])
        epoxide_carbons.add(match[2])
    
    # (5) Prepare helper functions for chain search.
    ring_info = mol.GetRingInfo()
    rings = [set(r) for r in ring_info.AtomRings()]
    
    def in_large_ring(atom_idx):
        # Return True if the atom is in any ring larger than 3 atoms.
        for r in rings:
            if atom_idx in r and len(r) > 3:
                return True
        return False
    
    def is_allowed_chain_atom(atom):
        # Only allow carbon atoms as part of the fatty acid chain.
        if atom.GetAtomicNum() != 6:
            return False
        idx = atom.GetIdx()
        # If not in any ring, it's allowed.
        if not atom.IsInRing():
            return True
        # If in a ring, allow only if every ring the atom is in is exactly 3 atoms.
        for r in rings:
            if idx in r and len(r) > 3:
                return False
        return True
    
    # (6) From an acid group, follow the chain.
    # We shall search from the acid carbon (first atom of acid match) along its allowed carbon neighbors.
    # Instead of simply measuring a maximum length, we will collect full paths.
    def dfs_paths(current_idx, visited, current_path):
        """
        A DFS that yields complete simple paths (lists of atom indices)
        following allowed chain atoms.
        """
        current_atom = mol.GetAtomWithIdx(current_idx)
        found_extension = False
        for nbr in current_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            if is_allowed_chain_atom(nbr):
                visited.add(nbr_idx)
                # Append neighbor index to the current path
                new_path = current_path + [nbr_idx]
                # Recurse further
                yield from dfs_paths(nbr_idx, visited, new_path)
                visited.remove(nbr_idx)
                found_extension = True
        if not found_extension:
            # If no extension, yield the current path.
            yield current_path

    # (7) For each acid match, try to find a chain path meeting our criteria.
    # We require: (a) at least 12 carbons in the chain (excluding the acid carbon),
    # and (b) one of the internal positions (not the first or last carbon of the path) is in epoxide_carbons.
    qualifying_path_found = False
    reason_details = []
    for match in acid_matches:
        acid_carbon_idx = match[0]  # the carbon of the acid carbonyl
        acid_atom = mol.GetAtomWithIdx(acid_carbon_idx)
        # Look at all allowed neighbors of the acid carbon.
        for nbr in acid_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and is_allowed_chain_atom(nbr):
                visited = {acid_carbon_idx, nbr.GetIdx()}
                # Start a DFS path with the first neighbor.
                for path in dfs_paths(nbr.GetIdx(), visited, [nbr.GetIdx()]):
                    # path is a list of carbon indices; we require a minimum chain length.
                    if len(path) < 12:
                        continue
                    # Check that the epoxide carbon is on this path (but not at the very beginning or very end).
                    # We only require one such occurrence.
                    internal_indices = path[1:-1]  # exclude first and last carbon in the path
                    if any(atom_idx in epoxide_carbons for atom_idx in internal_indices):
                        qualifying_path_found = True
                        break
                if qualifying_path_found:
                    break
        if qualifying_path_found:
            break

    if not qualifying_path_found:
        return False, ("Aliphatic chain too short or epoxide ring not "
                       "embedded in a long, linear chain starting from the acid group")
    
    return True, "Molecule has a carboxylic acid group, a sufficiently long linear aliphatic chain with an internal epoxide ring."

# For testing examples if this file is run as main.
if __name__ == "__main__":
    # Test one true positive example:
    test_smiles = "O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O"
    result, reason = is_epoxy_fatty_acid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)