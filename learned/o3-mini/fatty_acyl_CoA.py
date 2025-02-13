"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA
Definition:
  An acyl-CoA that results from the formal condensation of the thiol group 
  of coenzyme A with the carboxy group of any fatty acid.
  
Improved Strategy:
  1. Parse the SMILES.
  2. Look for a CoA fragment using a SMARTS pattern that targets a characteristic core 
     (avoiding the variable phosphate substituents).
  3. Look for a thioester bond (C(=O)S).
  4. For each thioester match, from the carbonyl carbon, choose a neighbor that is:
       • Not the carbonyl oxygen or bonded sulfur,
       • Not part of the CoA fragment,
       • And is a carbon atom.
  5. Use a recursive search (DFS) to traverse only contiguous aliphatic (non‐aromatic) 
     carbon atoms. This ensures that we only count the pure fatty acyl chain.
  6. Accept the molecule if one finds an acyl chain of at least four carbons.
     
This aims both to recognize valid fatty acyl-CoAs (including short-chain ones) and to 
reject acyl-CoAs whose “acyl” part contains extra heteroatoms or aromatic rings.
"""

from rdkit import Chem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if classified as fatty acyl-CoA, else False.
        str: Explanation for the classification.
    """
    # Parse the input
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a stricter SMARTS for a characteristic CoA core fragment.
    # This pattern is chosen to capture a part of the CoA moiety that is common while avoiding the variable phosphate groups.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C;R0]")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "Coenzyme A moiety not found"
    
    # Gather atom indices that are part of any CoA fragment
    coa_atom_indices = set()
    for match in coa_matches:
        coa_atom_indices.update(match)
    
    # Define SMARTS pattern for the thioester bond (carbonyl C(=O) bonded to an S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester bond (C(=O)S) not found; cannot be acyl-CoA"
    
    # Define a helper function (DFS) to compute the longest contiguous chain of aliphatic (non-aromatic) carbons.
    def longest_chain(atom, visited):
        # Only count carbon atoms that are not aromatic.
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            return 0
        max_length = 1  # count current carbon
        visited.add(atom.GetIdx())
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            # Skip if neighbor was visited or is part of the CoA core.
            if n_idx in visited or n_idx in coa_atom_indices:
                continue
            if neighbor.GetAtomicNum() == 6 and not neighbor.GetIsAromatic():
                branch_length = 1 + longest_chain(neighbor, visited.copy())
                if branch_length > max_length:
                    max_length = branch_length
        return max_length
    
    fatty_chain_found = False
    acyl_chain_length = 0
    
    # Iterate over each thioester fragment.
    # Each match returns a tuple (carbonyl, oxygen, sulfur).
    for match in thioester_matches:
        carbonyl_idx, oxygen_idx, sulfur_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # From the carbonyl, choose neighbor(s) that are not the oxygen, not the sulfur, and not part of CoA.
        acyl_neighbors = []
        for neighbor in carbonyl_atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in (oxygen_idx, sulfur_idx) or (n_idx in coa_atom_indices):
                continue
            if neighbor.GetAtomicNum() == 6:
                acyl_neighbors.append(neighbor)
        if not acyl_neighbors:
            continue  # try the next thioester match
        
        # For simplicity, use the first eligible neighbor as the start of the fatty acyl chain.
        start_atom = acyl_neighbors[0]
        chain_length = longest_chain(start_atom, set())
        # Accept if chain length is at least 4 carbon atoms.
        if chain_length >= 4:
            fatty_chain_found = True
            acyl_chain_length = chain_length
            break
    
    if not fatty_chain_found:
        return (False, "No sufficiently long fatty acyl chain (>=4 contiguous carbon atoms) detected")
    
    reason = (f"Contains a CoA moiety and a thioester bond linking an acyl chain "
              f"with {acyl_chain_length} contiguous carbon atoms (fatty acyl chain).")
    return True, reason