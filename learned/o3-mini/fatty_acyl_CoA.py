"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA
Definition:
  An acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any fatty acid.
  
The molecule must contain a CoA “tail” (we match a characteristic fragment found in most CoA derivatives)
and a thioester bond (C(=O)S) linking an acyl chain.
"""

from rdkit import Chem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    
    The strategy is:
      1. Parse the SMILES.
      2. Look for a CoA fragment. In many acyl-CoA compounds the fragment
         "SCCNC(=O)CCNC(=O)" is part of the CoA moiety.
      3. Look for a thioester bond ("C(=O)S"). Then find the acyl chain 
         attached to the carbonyl side of that group.
      4. Do a simple graph search (DFS) starting from the neighbor carbon (not the S)
         to count contiguous carbon atoms (ignoring atoms that are part of the CoA fragment).
         If the acyl chain has at least 3 carbons, we consider it a fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a fatty acyl-CoA, else False.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a characteristic CoA fragment.
    # This pattern is not exhaustive but should be found in many acyl-CoA compounds.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "Coenzyme A moiety not found"
    
    # Define a SMARTS pattern for a thioester bond (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester bond (C(=O)S) not found; cannot be acyl-CoA"
    
    # Get the set of atom indices that are part of the CoA fragment (from any match)
    coa_atom_indices = set()
    for match in coa_matches:
        for idx in match:
            coa_atom_indices.add(idx)
    
    # We write a simple DFS function that, starting from a given carbon atom,
    # walks through connected carbon atoms (ignoring any atoms that are in the CoA fragment)
    # and returns the number of carbons in that connected acyl chain.
    def dfs_count_carbons(atom, visited):
        count = 0
        # Only count if it is a carbon atom (atomic num 6)
        if atom.GetAtomicNum() == 6:
            count = 1
        else:
            return 0
        visited.add(atom.GetIdx())
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in visited or n_idx in coa_atom_indices:
                continue
            # Allow passage through carbon atoms (even if unsaturated)
            if neighbor.GetAtomicNum() == 6:
                count += dfs_count_carbons(neighbor, visited)
        return count

    # Now, among the thioester matches, we look for the acyl chain.
    # A thioester SMARTS match gives a tuple like (c_idx, o_idx, s_idx)
    fatty_chain_found = False
    acyl_chain_length = 0
    for match in thioester_matches:
        carbonyl_idx, oxygen_idx, sulfur_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # From the carbonyl atom, find the neighbor that is not the oxygen or sulfur.
        acyl_neighbors = []
        for neighbor in carbonyl_atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in (oxygen_idx, sulfur_idx) and n_idx not in coa_atom_indices:
                acyl_neighbors.append(neighbor)
        if not acyl_neighbors:
            continue  # This thioester match does not have an acyl chain (might be reversed)
        # For simplicity, take the first acyl neighbor and count contiguous carbons.
        start_atom = acyl_neighbors[0]
        count = dfs_count_carbons(start_atom, set())
        # We require the acyl chain to have at least 3 carbons.
        if count >= 3:
            fatty_chain_found = True
            acyl_chain_length = count
            break

    if not fatty_chain_found:
        return False, "No sufficiently long fatty acyl chain (>=3 contiguous carbons) detected"
    
    # Passed all checks, so classify as fatty acyl-CoA.
    reason = (f"Contains a CoA moiety, a thioester bond linking an acyl chain "
              f"with {acyl_chain_length} contiguous carbon(s) (fatty acyl chain).")
    return True, reason