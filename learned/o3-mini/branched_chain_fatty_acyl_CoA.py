"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any branched-chain fatty acid.
The improved algorithm:
  1. Checks for presence of a CoA moiety (via an adenine-derived substructure).
  2. Checks for a thioester bond linking the fatty acyl chain.
  3. From the thioester carbonyl, it identifies the fatty acyl neighbor (the carbon 
     not linked to S) and then performs a DFS that is restricted to carbon atoms,
     excluding the carbonyl atom, the sulfur attached to CoA, and any atoms that 
     belong to the known CoA substructure.
  4. The obtained acyl chain must be long enough (here, at least 3 carbons) and 
     at least one carbon in that set must have >2 carbon neighbors (within the set)
     indicating branching.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.

    The improved algorithm:
      1. Parses the SMILES.
      2. Requires an adenine substructure (used as a marker for CoA).
      3. Requires a thioester linkage [CX3](=O)[SX2].
      4. For each thioester match, identifies the fatty acyl chain by selecting the
         carbonyl's neighbor that is not sulfur (and not oxygen) and then performing a DFS
         restricted to carbon atoms. During the DFS, we do not traverse back to the carbonyl,
         the attached sulfur, or any atom that belongs to the CoA adenine substructure.
      5. Checks that the acyl chain is at least 3 carbons long and that at least one atom within
         the chain has more than 2 neighbors (from the chain) indicating branching.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a branched-chain fatty acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Step 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Step 2. Check for the presence of a CoA moiety by matching the adenine substructure.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"
    # Get indices for all adenine matches (we combine all atoms that were found)
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    coa_atoms = set()
    for match in coa_matches:
        coa_atoms.update(match)
    
    # Step 3. Look for the thioester linkage using SMARTS pattern.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester bond not found"

    branched_found = False
    linear_found = False  # At least one valid fatty acyl chain was found but appears linear
    # Process each thioester match
    for match in thioester_matches:
        # match[0] is the carbonyl carbon; match[1] is the sulfur atom attached to CoA.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Step 4. Identify the fatty acyl neighbor: a neighbor of the carbonyl that is not the sulfur
        # and not the double-bonded oxygen.
        fatty_acyl_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue  # skip the sulfur atom that goes to CoA
            if nbr.GetAtomicNum() == 8:  # skip oxygen (from =O)
                continue
            # We assume this neighbor initiates the fatty acyl chain.
            fatty_acyl_neighbor = nbr
            break
        if fatty_acyl_neighbor is None:
            continue  # try next thioester match if none found

        # Step 5. Perform a DFS starting at fatty_acyl_neighbor gathering only carbon atoms.
        # Exclude the carbonyl atom, the sulfur, and any atoms belonging to CoA.
        visited = set()
        stack = [fatty_acyl_neighbor.GetIdx()]
        while stack:
            current_idx = stack.pop()
            if current_idx in visited:
                continue
            # Do not traverse into the carbonyl or the sulfur.
            if current_idx == carbonyl_idx or current_idx == sulfur_idx:
                continue
            # Also do not go into atoms that are part of the CoA moiety.
            if current_idx in coa_atoms:
                continue
            atom = mol.GetAtomWithIdx(current_idx)
            if atom.GetAtomicNum() != 6:
                continue  # only consider carbon atoms
            visited.add(current_idx)
            # Traverse all neighbors that are carbons.
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in visited:
                    continue
                # Avoid going back to the carbonyl or sulfur or into the CoA region.
                if nbr_idx == carbonyl_idx or nbr_idx == sulfur_idx or nbr_idx in coa_atoms:
                    continue
                if nbr.GetAtomicNum() == 6:
                    stack.append(nbr_idx)
        acyl_chain = visited

        # Check the chain length (at least 3 carbons).
        if len(acyl_chain) < 3:
            continue  # chain too short; try next thioester match

        # Step 6. Check for branching: in a perfectly linear chain,
        # every internal carbon (nonterminal) has exactly 2 carbon neighbors (within the acyl chain)
        # and terminals have one. Therefore, a branch is present if any carbon atom in the chain
        # has more than 2 neighbors (from the chain).
        branch_found = False
        for idx in acyl_chain:
            atom = mol.GetAtomWithIdx(idx)
            # Count neighbors that are in the acyl chain (ignore others).
            n_chain_neighbors = sum(1 for nbr in atom.GetNeighbors() if nbr.GetIdx() in acyl_chain)
            if n_chain_neighbors > 2:
                branch_found = True
                break

        if branch_found:
            branched_found = True
            # If one match gives a branch, we conclude it is a branched-chain fatty acyl-CoA.
            return True, "Branched-chain fatty acyl-CoA detected"
        else:
            linear_found = True
            # Do not return immediately; try other matches.
    
    # End of loop processing all thioester matches.
    if linear_found:
        return False, "Acyl chain appears linear; no branch detected"
    else:
        return False, "No fatty acyl chain meeting criteria found"


# Example usage:
if __name__ == "__main__":
    # Test with one of the provided SMILES strings.
    sample_smiles = "CC(C)CCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_branched_chain_fatty_acyl_CoA(sample_smiles)
    print(result, reason)