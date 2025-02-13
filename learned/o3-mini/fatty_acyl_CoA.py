"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA 
Definition:
  An acyl-CoA that results from the formal condensation of the thiol group 
  of coenzyme A with the carboxy group of any fatty acid.

The molecule must contain a CoA moiety (we search for a characteristic fragment)
and a thioester bond (C(=O)S) linking an acyl chain. In our improved classifier we:
  • Use a stricter SMARTS to find a core CoA fragment.
  • Identify the thioester bond.
  • From the carbonyl carbon, choose a neighbor that is not the oxygen or the sulfur 
    of the thioester and not in the CoA fragment.
  • Use a recursive chain search that only accepts contiguous aliphatic carbon atoms (atomic number 6).
  • Require that the longest chain is at least 5 atoms long.
"""

from rdkit import Chem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    
    Improved Strategy:
      1. Parse the SMILES.
      2. Look for a CoA fragment using a SMARTS that matches a characteristic core 
         (we try to avoid variable phosphate substituents).
      3. Look for a thioester bond (C(=O)S).
      4. For each thioester, examine the carbonyl carbon’s neighbors. Exclude the oxygen 
         (of the C=O), the connecting sulfur, and any atom in the CoA fragment.
         Also require that the chosen neighbor is a carbon.
      5. Recursively compute the longest contiguous chain of aliphatic carbons (atomic num 6)
         allowing only these atoms. (No sulfur allowed in the acyl chain now.)
      6. Accept if a chain of at least 5 carbons is found.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a fatty acyl-CoA, else False.
        str: Explanation for the classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use a stricter SMARTS for a characteristic CoA core fragment.
    # (This pattern covers a fragment that is common in many acyl-CoA compounds.)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C;R0]")  
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "Coenzyme A moiety not found"
    
    # Collect atom indices for the CoA fragment(s)
    coa_atom_indices = set()
    for match in coa_matches:
        coa_atom_indices.update(match)
    
    # SMARTS pattern for a thioester bond: carbonyl (C(=O)) bonded to an S.
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester bond (C(=O)S) not found; cannot be acyl-CoA"
    
    # Define a DFS to compute the longest contiguous chain consisting only of aliphatic carbons.
    # We only allow carbon (atomic number 6) in the fatty acyl chain.
    def longest_chain(atom, visited):
        if atom.GetAtomicNum() != 6:
            return 0
        max_length = 1  # count this carbon 
        visited.add(atom.GetIdx())
        # Explore neighbors that are carbon and not already visited and not in the CoA fragment.
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in visited or n_idx in coa_atom_indices:
                continue
            if neighbor.GetAtomicNum() == 6:
                branch_length = 1 + longest_chain(neighbor, visited.copy())
                if branch_length > max_length:
                    max_length = branch_length
        return max_length

    fatty_chain_found = False
    acyl_chain_length = 0

    # Process each thioester match.
    # Each match is a tuple of indices: (carbonyl, oxygen, sulfur).
    for match in thioester_matches:
        carbonyl_idx, oxygen_idx, sulfur_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # From carbonyl atom, look for neighbors that are neither the oxygen nor the attached sulfur,
        # and not part of the CoA fragment.
        acyl_neighbors = []
        for neighbor in carbonyl_atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in (oxygen_idx, sulfur_idx) or n_idx in coa_atom_indices:
                continue
            # Require that the neighbor is a carbon (to start a fatty acyl chain)
            if neighbor.GetAtomicNum() == 6:
                acyl_neighbors.append(neighbor)
        if not acyl_neighbors:
            continue  # No suitable acyl chain found at this thioester bond

        # For simplicity, choose the first candidate as the acyl chain start.
        start_atom = acyl_neighbors[0]
        chain_length = longest_chain(start_atom, set())
        if chain_length >= 5:
            fatty_chain_found = True
            acyl_chain_length = chain_length
            break
    
    if not fatty_chain_found:
        return (False, 
                "No sufficiently long fatty acyl chain (>=5 contiguous carbon atoms) detected")
    
    reason = (f"Contains a CoA moiety, a thioester bond linking an acyl chain "
              f"with {acyl_chain_length} contiguous carbon atoms (fatty acyl chain).")
    return True, reason