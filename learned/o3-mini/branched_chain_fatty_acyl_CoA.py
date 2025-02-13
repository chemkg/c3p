"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any branched-chain fatty acid.
The improved algorithm:
  1. Checks for a CoA marker (via an adenine-derived substructure).
  2. Locates a thioester linkage [CX3](=O)[SX2].
  3. From the thioester carbonyl, it identifies the fatty acyl neighbor (the atom not linked to S and not oxygen)
     and then collects the contiguous set of carbon atoms (via DFS) that are not part of CoA, the carbonyl, or the sulfur.
  4. The obtained acyl chain must be at least 3 carbons long and at least one carbon in that set must have at least 3 neighbors 
     (counting only neighbors in the chain) to indicate branching.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.

    The improved algorithm:
      1. Parses the molecule.
      2. Checks for the presence of a CoA moiety (by matching an adenine substructure).
      3. Finds a thioester bond ([CX3](=O)[SX2]) bridging the fatty acyl chain with the CoA.
      4. For each thioester match, identifies the fatty acyl's starting atom—the carbon adjacent to the carbonyl
         (which is not the sulfur and not an oxygen).
      5. Using a depth-first search (DFS), collects all connected carbon atoms (atomic number 6) that belong to
         the acyl chain while forbidding reentry into the carbonyl, the sulfur, or any atoms already in the CoA moiety.
      6. If the resulting set has at least 3 carbon atoms and if any carbon in that subset has three or more neighbors
         (when only neighbors in the subset are counted), then a branch is present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is identified as a branched-chain fatty acyl-CoA.
        str: Explanation for the classification.
    """
    # Step 1. Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2. Check for a CoA moiety (use adenine SMARTS as a marker)
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"
    # Get all atoms that match the CoA adenine ring
    coa_atoms = set()
    for match in mol.GetSubstructMatches(coa_pattern):
        coa_atoms.update(match)
    
    # Step 3. Locate thioester linkage [CX3](=O)[SX2]
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester bond not found"
    
    branched_found = False
    linear_found = False  # at least one fatty acyl chain was found but appears linear

    # Process each thioester match
    for match in thioester_matches:
        # match[0] is the carbonyl carbon; match[1] is the sulfur directly attached
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Step 4. Identify the fatty acyl neighbor from the carbonyl:
        # (choose a neighbor that is not the sulfur atom and not an oxygen [=O])
        fatty_acyl_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue
            if nbr.GetAtomicNum() == 8:  # skip oxygen atom(s)
                continue
            fatty_acyl_neighbor = nbr
            break
        if fatty_acyl_neighbor is None:
            continue  # no valid fatty acyl neighbor in this thioester, try next
        
        # Step 5. From fatty_acyl_neighbor, perform a DFS collecting only carbon atoms that belong to the acyl chain.
        # We exclude atoms that are the carbonyl, the sulfur or part of the CoA (coa_atoms).
        chain_set = set()
        stack = [fatty_acyl_neighbor.GetIdx()]
        while stack:
            current_idx = stack.pop()
            if current_idx in chain_set:
                continue
            if current_idx == carbonyl_idx or current_idx == sulfur_idx:
                continue
            if current_idx in coa_atoms:
                continue
            atom = mol.GetAtomWithIdx(current_idx)
            if atom.GetAtomicNum() != 6:
                continue  # only accept carbon atoms
            chain_set.add(current_idx)
            # Now add neighbors that are carbons and not disallowed.
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in chain_set:
                    continue
                if nbr_idx == carbonyl_idx or nbr_idx == sulfur_idx or nbr_idx in coa_atoms:
                    continue
                if nbr.GetAtomicNum() == 6:
                    stack.append(nbr_idx)
        
        # Check if the acyl chain is long enough
        if len(chain_set) < 3:
            continue  # too short, try another thioester match
        
        # Step 6. Check for branching.
        # Here we count, for each atom in the chain_set, the number of neighbors (also in chain_set).
        branch_here = False
        for idx in chain_set:
            atom = mol.GetAtomWithIdx(idx)
            chain_neighbors = 0
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in chain_set:
                    chain_neighbors += 1
            # In a linear chain, internal carbons have 2 neighbors and terminal carbons have 1.
            # A branch point should yield at least 3 neighbors.
            if chain_neighbors > 2:
                branch_here = True
                break
        
        if branch_here:
            branched_found = True
            return True, "Branched-chain fatty acyl-CoA detected"
        else:
            linear_found = True  # valid acyl chain was found but no branch
    
    # End of processing thioester matches.
    if linear_found:
        return False, "Acyl chain appears linear; no branch detected"
    else:
        return False, "No fatty acyl chain meeting criteria found"


# Example usage (for testing):
if __name__ == "__main__":
    # Example SMILES from the provided list: isoheptadecanoyl-CoA
    test_smiles = "CC(C)CCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_branched_chain_fatty_acyl_CoA(test_smiles)
    print(result, reason)