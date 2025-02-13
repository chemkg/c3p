"""
Classifies: CHEBI:25029 leukotriene
"""
#!/usr/bin/env python
"""
Classifies: leukotriene
Definition:
  A leukotriene is any icosanoid derived from arachidonic acid (a C20 polyunsaturated fatty acid)
  produced by leukocytes that features a backbone of 20 carbons containing exactly 4 carbon–carbon 
  double bonds, and at least one contiguous segment showing a D–S–D–S–D pattern (i.e. alternating 
  double and single bonds over 5 bonds).
  
Our approach is to:
  1. Parse the molecule.
  2. Identify “allowed” carbon atoms – those that are not in rings of size >3 (we allow 3‐membered rings 
     to capture epoxide features).
  3. Build a connectivity (neighbors) dictionary among allowed carbons.
  4. Use DFS to find the longest simple (acyclic) chain through the allowed carbons.
  5. Slide a window of exactly 20 carbons along this chain. For each candidate chain, examine the bonds 
     between successive carbons. Accept the candidate if it has exactly 4 C=C bonds and at least one 5‐bond 
     segment with the pattern D–S–D–S–D.
     
If any candidate 20‐carbon backbone is found that meets the criteria, the molecule is classified as a leukotriene.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    
    Criteria:
      1. The molecule must be parsable.
      2. It must contain a contiguous subchain of exactly 20 acyclic carbon atoms (ignoring carbons in rings >3)
         that forms the “backbone.”
      3. Along the backbone there must be exactly 4 carbon–carbon double bonds.
      4. Within those bonds there must be at least one contiguous segment with the pattern D–S–D–S–D,
         where D indicates a double bond and S a single bond.
    
    Args:
      smiles (str): SMILES string for the molecule.
    
    Returns:
      (bool, str): Tuple with (True, reason) if classified as a leukotriene;
                   (False, reason) if not; may return (None, None) if task is infeasible.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain ring information.
    ring_info = mol.GetRingInfo().AtomRings()
    # Create a set of atom indices that belong to rings of size > 3.
    non_epoxide_ring_atoms = set()
    for ring in ring_info:
        if len(ring) > 3:
            non_epoxide_ring_atoms.update(ring)
    
    # Allowed carbons: only carbon atoms (atomic number 6) that are not in rings of size > 3.
    allowed_carbon_idxs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            if atom.GetIdx() in non_epoxide_ring_atoms:
                continue
            allowed_carbon_idxs.append(atom.GetIdx())
    
    if not allowed_carbon_idxs:
        return False, "No suitable acyclic carbon atoms found (or all carbons are in rings)"
    
    # Build neighbor dictionary for allowed carbons.
    neighbors = {idx: [] for idx in allowed_carbon_idxs}
    for idx in allowed_carbon_idxs:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in allowed_carbon_idxs:
                neighbors[idx].append(nbr.GetIdx())
    
    # Use DFS to find one of the longest simple paths (i.e., no repeated atoms)
    longest_path = []
    total_allowed = len(allowed_carbon_idxs)
    
    def dfs(current, path, visited):
        nonlocal longest_path
        # Bound: if maximum possible chain length from here cannot exceed current longest, prune search.
        if len(path) + (total_allowed - len(visited)) <= len(longest_path):
            return
        if len(path) > len(longest_path):
            longest_path = path[:]
        for nbr in neighbors.get(current, []):
            if nbr not in visited:
                visited.add(nbr)
                dfs(nbr, path + [nbr], visited)
                visited.remove(nbr)
    
    for start in allowed_carbon_idxs:
        dfs(start, [start], set([start]))
        if len(longest_path) == total_allowed:
            break  # maximum possible chain found
    
    # If the longest chain is shorter than 20, then it cannot be a leukotriene.
    if len(longest_path) < 20:
        return False, f"Longest acyclic carbon chain has {len(longest_path)} carbons; expected at least 20."
    
    # Now, slide a window of 20 carbons along the longest chain to find a candidate backbone.
    def check_candidate(candidate):
        # candidate is a list of 20 atom indices.
        bond_types = []  # record "D" for double bonds, "S" for single bonds.
        double_bond_count = 0
        for i in range(len(candidate) - 1):
            bond = mol.GetBondBetweenAtoms(candidate[i], candidate[i+1])
            if bond is None:
                bond_types.append("NA")
                continue
            # If the bond is a double bond we mark it as "D".
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                bond_types.append("D")
                double_bond_count += 1
            else:
                bond_types.append("S")
        if double_bond_count != 4:
            return False, f"Candidate backbone has {double_bond_count} C=C bonds; expected exactly 4."
        
        # Look for at least one contiguous window of 5 bonds that matches D-S-D-S-D.
        found_conjugation = False
        for i in range(len(bond_types) - 4):
            segment = bond_types[i:i+5]
            if segment == ["D", "S", "D", "S", "D"]:
                found_conjugation = True
                break
        if not found_conjugation:
            return False, "No contiguous conjugated tri-double-bond system (D-S-D-S-D) found in candidate backbone."
        return True, f"Candidate backbone (20 carbons) has exactly 4 C=C bonds and a D-S-D-S-D conjugated segment."
    
    # Try every contiguous subchain of length 20 in the longest chain.
    for i in range(0, len(longest_path) - 20 + 1):
        candidate = longest_path[i:i+20]
        valid, reason = check_candidate(candidate)
        if valid:
            return True, (f"Found a 20-carbon backbone (from atom indices {candidate[0]} to {candidate[-1]}) with exactly "
                          f"4 C=C bonds and a conjugated segment. {reason}")
    
    # If none of the contiguous subchains qualifies, report failure.
    return False, "No 20-carbon subchain with exactly 4 C=C bonds and a D-S-D-S-D conjugated segment was found."

# Example usage (uncomment to test):
# test_smiles = "CCCC\\C=C/C[C@@H](O)\\C=C\\C=C\\C=C\\[C@@H](O)CCCC(O)=O"  # 6-trans-leukotriene B4
# print(is_leukotriene(test_smiles))