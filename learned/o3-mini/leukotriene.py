"""
Classifies: CHEBI:25029 leukotriene
"""
#!/usr/bin/env python
"""
Classifies: leukotriene
Definition:
  A leukotriene is defined as any icosanoid derived from arachidonic acid (a C20 polyunsaturated fatty acid)
  and its derivatives that features a backbone with 20 carbons containing exactly 4 carbon–carbon double bonds,
  and in which at least one segment of three alternating (conjugated) double bonds (D–S–D–S–D) is present.
  
We attempt to identify the “backbone” as the longest continuous (acyclic) chain among carbons that are not in rings,
except that carbons in small (3–membered) rings (which may be present in epoxides) are allowed.
"""

from rdkit import Chem

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    
    The molecule must satisfy:
      1. It must be parseable by RDKit.
      2. Its longest continuous carbon backbone (ignoring carbons in rings of size >3) must have at least 20 carbons.
      3. There must be exactly 4 carbon–carbon double bonds along that backbone.
      4. There must be at least one contiguous segment showing a D–S–D–S–D bond pattern along that backbone,
         where 'D' means a C=C double bond and 'S' a single bond.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): (True, reason) if the molecule meets the criteria, (False, reason) otherwise.
                  If the task times out or cannot be done, (None, None) may be returned.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information.
    ring_info = mol.GetRingInfo().AtomRings()
    # Build a set of atom indices that are in any ring with size > 3.
    non_epoxide_ring_atoms = set()
    for ring in ring_info:
        if len(ring) > 3:
            for idx in ring:
                non_epoxide_ring_atoms.add(idx)
    
    # We now want to consider only carbon atoms (atomic number 6) that are not in a non-epoxide ring.
    allowed_carbon_idxs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            # If the carbon is in a non-3-membered ring, skip it.
            if atom.GetIdx() in non_epoxide_ring_atoms:
                continue
            allowed_carbon_idxs.append(atom.GetIdx())
    
    if not allowed_carbon_idxs:
        return False, "No suitable acyclic carbon atoms found (or all carbons are in rings)"
    
    # Build a neighbor dictionary for allowed carbon atoms.
    neighbors = {idx: [] for idx in allowed_carbon_idxs}
    for idx in allowed_carbon_idxs:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in allowed_carbon_idxs:
                # Only add if the neighbor is in our allowed set.
                neighbors[idx].append(nbr.GetIdx())
    
    # Use DFS with branch-and-bound to find the longest simple path (only through allowed carbons).
    longest_path = []
    total_allowed = len(allowed_carbon_idxs)
    
    def dfs(current, path, visited):
        nonlocal longest_path
        # Bound: maximum possible length from here cannot exceed current_length + remaining nodes.
        if len(path) + (total_allowed - len(visited)) <= len(longest_path):
            return
        # Update longest path if needed.
        if len(path) > len(longest_path):
            longest_path = path[:]
        for nbr in neighbors.get(current, []):
            if nbr not in visited:
                visited.add(nbr)
                dfs(nbr, path + [nbr], visited)
                visited.remove(nbr)
    
    # Run DFS starting from each allowed carbon.
    for start in allowed_carbon_idxs:
        dfs(start, [start], set([start]))
        # If the longest possible chain among allowed carbons equals total allowed, no need to search further.
        if len(longest_path) == total_allowed:
            break

    if len(longest_path) < 20:
        return False, f"Longest carbon chain in acyclic region has {len(longest_path)} carbons; expected at least 20."

    # Now, analyze the bonds along the found chain.
    # We record bond types between consecutive carbons in the chain.
    chain_bond_types = []
    double_bond_count = 0
    for i in range(len(longest_path)-1):
        bond = mol.GetBondBetweenAtoms(longest_path[i], longest_path[i+1])
        if bond is None:
            chain_bond_types.append("NA")
            continue
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count += 1
            chain_bond_types.append("D")
        else:
            chain_bond_types.append("S")
    
    if double_bond_count != 4:
        return False, f"Found {double_bond_count} C=C bonds along the backbone; expected exactly 4."
    
    # Check for at least one contiguous D-S-D-S-D segment.
    found_conjugated = False
    # We need a sliding window of 5 bond types. 
    for i in range(len(chain_bond_types) - 4):
        segment = chain_bond_types[i:i+5]
        if segment == ["D", "S", "D", "S", "D"]:
            found_conjugated = True
            break
    if not found_conjugated:
        return False, "No contiguous conjugated tri-double-bond system (D-S-D-S-D) found along the backbone."
    
    return True, (f"Molecule has a backbone of {len(longest_path)} acyclic carbons, "
                  f"with exactly 4 C=C bonds and a conjugated (D-S-D-S-D) segment typical of leukotrienes.")

# Example test (uncomment to run a test case):
# test_smiles = "CCCC\\C=C/C[C@@H](O)\\C=C\\C=C\\C=C\\[C@@H](O)CCCC(O)=O"  # 6-trans-leukotriene B4 example
# print(is_leukotriene(test_smiles))