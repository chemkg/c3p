"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any branched-chain fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    The method uses the following heuristics:
      1. It must contain a CoA portion (detected by an adenine substructure).
      2. It must contain a thioester group (C(=O)S) linking the acyl moiety to CoA.
      3. The fatty acyl chain (the part attached to the carbonyl carbon not going to the S)
         is traversed (limiting to carbon atoms) and checked for branching (any carbon
         with more than two connections within the chain indicates branching).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a branched-chain fatty acyl-CoA, False otherwise.
        str: Reason for the classification.
    """
    
    # Step 1. Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2. Check for the presence of a CoA moiety.
    # Here we look for an adenine substructure commonly found in CoA.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"
    
    # Step 3. Check for the thioester group linking fatty acyl and CoA.
    # This SMARTS matches a carbonyl carbon attached to a sulfur: C(=O)S.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester bond not found"
    
    # We loop through each thioester match to see if any gives a branched-chain acyl part.
    for match in thioester_matches:
        # In the pattern match, match[0] is the carbonyl C and match[1] is the S.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the neighbor of the carbonyl carbon that is not the thioester S 
        # nor the oxygen from the carbonyl (which is implicit/double-bonded).
        fatty_acyl_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue  # skip the sulfur leading to the CoA portion
            # Skip oxygen (likely the double-bonded oxygen in the carbonyl)
            if nbr.GetAtomicNum() == 8:
                continue
            # Assume that this neighbor (usually carbon) is the beginning of the fatty acyl chain.
            fatty_acyl_neighbor = nbr
            break
        if fatty_acyl_neighbor is None:
            continue  # try the next thioester match
        
        # Step 4. Traverse the acyl chain (starting from the carbonyl carbon) restricting to carbon atoms.
        # We exclude atoms from the CoA portion (the sulfur branch) and any heteroatoms.
        acyl_atoms = set()
        stack = [carbonyl_atom]
        while stack:
            a = stack.pop()
            # Only follow carbon atoms (atomic num 6) to build the fatty acyl chain.
            if a.GetIdx() in acyl_atoms:
                continue
            if a.GetAtomicNum() != 6:
                continue
            acyl_atoms.add(a.GetIdx())
            for nbr in a.GetNeighbors():
                # Avoid leaving the acyl chain by not following bonds to non-carbons.
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in acyl_atoms:
                    # Also (if desired) one might check bond type for further restrictions.
                    stack.append(nbr)
        
        # Require a minimum chain length (for a fatty acyl chain)
        if len(acyl_atoms) < 3:
            continue  # try next match
        
        # Step 5. Determine if the acyl chain is branched.
        # In a linear alkane chain, internal carbons (not at the ends) have exactly 2 carbon neighbors,
        # while terminal carbons have 1. If any carbon in the acyl chain has >2 neighbors (within the chain),
        # then this chain is branched.
        branch_found = False
        for idx in acyl_atoms:
            atom = mol.GetAtomWithIdx(idx)
            # Count neighbors that are part of the acyl chain.
            n_acyl_neighbors = sum(1 for nbr in atom.GetNeighbors() if nbr.GetIdx() in acyl_atoms)
            if n_acyl_neighbors > 2:
                branch_found = True
                break
        
        if not branch_found:
            # Even if a thioester is present, if no branching is found then this is a linear fatty acyl.
            return False, "Acyl chain appears linear; no branch detected"
        else:
            return True, "Branched-chain fatty acyl-CoA detected"
    
    return False, "No fatty acyl chain meeting criteria found"
    
# Example usage:
if __name__ == "__main__":
    # You can test with one of the sample SMILES strings.
    sample_smiles = "CC(C)CCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_branched_chain_fatty_acyl_CoA(sample_smiles)
    print(result, reason)