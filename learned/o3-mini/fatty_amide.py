"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: Fatty Amide
Definition: A monocarboxylic acid amide derived from a fatty acid.
A fatty amide is recognized by an amide bond (C(=O)N) where the acyl portion (the carbonyl side)
comes from a fatty acid (i.e. contains a long aliphatic chain).
The heuristic used here is that the acyl chain must contain at least 8 carbons (including the carbonyl carbon).
"""

from rdkit import Chem

def longest_aliphatic_chain_length(atom, visited):
    """
    Recursively computes the longest length (in number of carbon atoms) 
    starting from the given atom. Traverses only through non‐aromatic, non‐ring carbon atoms.
    visited is a set of visited atom indices to avoid cycles.
    """
    max_length = 1  # count current atom
    # We use a copy of visited for different DFS branches.
    for nbr in atom.GetNeighbors():
        if (nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic() and not nbr.IsInRing() 
            and nbr.GetIdx() not in visited):
            new_visited = visited.copy()
            new_visited.add(nbr.GetIdx())
            length = 1 + longest_aliphatic_chain_length(nbr, new_visited)
            if length > max_length:
                max_length = length
    return max_length

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a fatty amide.
    
    A fatty amide is defined as a monocarboxylic acid amide derived from a fatty acid,
    meaning that there is an amide group (C(=O)N) such that the acyl (CO–) part
    is derived from a fatty acid (i.e. has a long aliphatic chain).
    
    For our heuristic, we require that the acyl chain (starting from the carbonyl carbon and
    following carbon-only, non-aromatic, acyclic paths) has at least 8 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the fatty amide criteria, False otherwise.
        str: A message giving the reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the amide substructure: C(=O)N
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    if amide_smarts is None:
        return False, "Error generating amide SMARTS pattern"

    # Find all substructure matches for the amide
    matches = mol.GetSubstructMatches(amide_smarts)
    if not matches:
        return False, "No amide group found in the molecule"

    # For each amide match, check whether the acyl part (the fatty acid–derived side) is long enough.
    # The SMARTS "C(=O)N" returns a tuple of atom indices (carbonyl carbon, oxygen, and nitrogen).
    for match in matches:
        # Extract atoms using the indices from the match.
        carbonyl = mol.GetAtomWithIdx(match[0])
        oxy = mol.GetAtomWithIdx(match[1])
        amide_nitrogen = mol.GetAtomWithIdx(match[2])
        
        # From the carbonyl carbon, get the neighbor(s) that are part of the acyl chain.
        # Exclude the oxygen (double-bonded) and the amide nitrogen.
        acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors() 
                          if nbr.GetIdx() not in {oxy.GetIdx(), amide_nitrogen.GetIdx()}
                          and nbr.GetAtomicNum() == 6
                          and not nbr.GetIsAromatic() 
                          and not nbr.IsInRing()]
        if not acyl_neighbors:
            # No carbon neighbors representing a fatty acyl chain found; check next match.
            continue

        # For each possible acyl connection, compute longest chain length.
        # We include the carbonyl carbon itself (chain length = 1 + chain from neighbor).
        max_chain_length = 0
        for nbr in acyl_neighbors:
            # Start DFS from the neighbor, making sure we do not go back to the carbonyl carbon.
            chain_length = 1 + longest_aliphatic_chain_length(nbr, visited={carbonyl.GetIdx()})
            if chain_length > max_chain_length:
                max_chain_length = chain_length
        
        # Use a threshold for the fatty acyl chain; here we require at least 8 carbon atoms.
        if max_chain_length >= 8:
            return True, f"Found fatty amide: amide group with acyl chain length {max_chain_length} carbons."
    
    return False, "Amide group found but no acyl chain long enough to be fatty (threshold 8 carbons)."