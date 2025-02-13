"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA 
Definition:
  An acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any fatty acid.
  
The molecule must contain a CoA “tail” (we match a characteristic fragment found in many CoA derivatives)
and a thioester bond (C(=O)S) linking an acyl chain.
We now improve the search for the acyl chain by – 
  • computing the longest contiguous chain of “aliphatic” atoms (carbon or sulfur)
  • requiring a minimum chain length of 5.
"""

from rdkit import Chem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    
    Strategy:
      1. Parse the SMILES.
      2. Look for a CoA fragment. Here we use a characteristic fragment that is found in many acyl-CoA compounds.
      3. Look for a thioester bond (C(=O)S).
      4. For a given thioester bond, from the carbonyl carbon, select the neighbor that is not the oxygen (or the sulfur of the thioester) nor part of the CoA fragment.
         Then use a recursive DFS to compute the maximum length of contiguous atoms in the acyl chain. Here we recursively allow only aliphatic atoms:
         atoms that are either carbon (atomic number 6) or sulfur (atomic number 16). 
      5. If the chain length is at least 5, we classify the molecule as a fatty acyl-CoA.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a fatty acyl-CoA, else False.
        str: Explanation for the classification result.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a characteristic CoA fragment.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "Coenzyme A moiety not found"
    
    # Collect all atom indices that are part of a CoA fragment (from any match).
    coa_atom_indices = set()
    for match in coa_matches:
        for idx in match:
            coa_atom_indices.add(idx)
    
    # SMARTS pattern for thioester bond: C(=O)S 
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester bond (C(=O)S) not found; cannot be acyl-CoA"
    
    # Define a recursive function to compute the longest contiguous acyl chain length.
    # We allow atoms that are carbon (atomic num 6) or sulfur (atomic num 16).
    # Instead of summing counts over branches we return the maximum linear chain length.
    def longest_chain(atom, visited):
        # Start with length 1 for this atom (if allowed)
        if atom.GetAtomicNum() not in (6, 16):
            return 0
        max_length = 1
        visited.add(atom.GetIdx())
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in visited or n_idx in coa_atom_indices:
                continue
            if neighbor.GetAtomicNum() in (6, 16):
                # Copy visited for each branch so as not to mix paths.
                branch_length = 1 + longest_chain(neighbor, visited.copy())
                if branch_length > max_length:
                    max_length = branch_length
        return max_length

    # Check each thioester match to try and locate an acyl chain.
    # thioester_matches gives tuples (carbonyl, oxygen, sulfur)
    fatty_chain_found = False
    acyl_chain_length = 0
    for match in thioester_matches:
        carbonyl_idx, oxygen_idx, sulfur_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Look at neighbors of the carbonyl atom that are not the oxygen (of C=O) and not the sulfur of the thioester,
        # and also not part of the CoA fragment.
        acyl_neighbors = []
        for neighbor in carbonyl_atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in (oxygen_idx, sulfur_idx) or n_idx in coa_atom_indices:
                continue
            acyl_neighbors.append(neighbor)
        if not acyl_neighbors:
            continue  # This thioester does not expose an acyl chain.
        # For simplicity, choose the first neighbor as the start of the acyl chain.
        start_atom = acyl_neighbors[0]
        chain_length = longest_chain(start_atom, set())
        # Require a minimum chain length of 5 aliphatic atoms.
        if chain_length >= 5:
            fatty_chain_found = True
            acyl_chain_length = chain_length
            break

    if not fatty_chain_found:
        return False, "No sufficiently long fatty acyl chain (>=5 contiguous aliphatic atoms) detected"
    
    reason = (f"Contains a CoA moiety, a thioester bond linking an acyl chain "
              f"with {acyl_chain_length} contiguous aliphatic atoms (fatty acyl chain).")
    return True, reason