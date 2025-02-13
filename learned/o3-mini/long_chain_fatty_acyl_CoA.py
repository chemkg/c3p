"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: long-chain fatty acyl-CoA
(Defined as a fatty acyl-CoA resulting from condensation of a fatty acid (with 13-22 carbons, including the carbonyl)
 with coenzyme A.)
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    
    The algorithm checks for:
      1. A CoA moiety by detecting a characteristic pantetheine fragment.
      2. A thioester linkage group ([CX3](=O)[SX2]) connecting the fatty acyl chain to CoA.
      3. A fatty acyl chain (starting from the carbonyl carbon) that contains between 13 and 22 carbon atoms.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a long-chain fatty acyl-CoA, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Detect the CoA moiety.
    # We check for the pantetheine fragment "SCCNC(=O)CCNC(=O)" that is characteristic of CoA.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_frag = Chem.MolFromSmarts(coa_smarts)
    if coa_frag is None:
        return False, "Invalid CoA SMARTS pattern"
    if not mol.HasSubstructMatch(coa_frag):
        return False, "CoA moiety not detected (missing pantetheine fragment)"
    
    # 2. Detect the thioester group linking the fatty acyl chain to CoA.
    # The SMARTS "[CX3](=O)[SX2]" should match three atoms: the carbonyl carbon, its double-bonded oxygen and the sulfur.
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if thioester is None:
        return False, "Invalid thioester SMARTS pattern"
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester group linking fatty acyl chain to CoA detected"
    
    reasons = []
    # 3. For each thioester group, try to isolate the fatty acyl chain attached.
    for match in thioester_matches:
        # Expect the match to contain three atoms: carbonyl carbon, oxygen, and sulfur.
        if len(match) != 3:
            reasons.append("Encountered a thioester match with unexpected atom count")
            continue
        carbonyl_idx, oxygen_idx, sulfur_idx = match  # Unpack the three indices.
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the neighbor of the carbonyl that is not the sulfur.
        # This neighbor is assumed to be the beginning of the fatty acyl chain.
        acyl_neighbors = []
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetIdx() == sulfur_idx:
                continue
            # We expect the acyl chain to consist mostly of carbon atoms.
            if neighbor.GetAtomicNum() == 6:
                acyl_neighbors.append(neighbor)
        if len(acyl_neighbors) != 1:
            reasons.append("Ambiguous or missing acyl chain for one thioester match")
            continue
        acyl_start_idx = acyl_neighbors[0].GetIdx()
        
        # Perform a depth-first search (DFS) along carbon–carbon bonds to count the number of carbons in the chain.
        # We include the carbonyl carbon in the count.
        acyl_atoms = set()
        acyl_atoms.add(carbonyl_idx)
        stack = [acyl_start_idx]
        while stack:
            curr_idx = stack.pop()
            if curr_idx in acyl_atoms:
                continue
            curr_atom = mol.GetAtomWithIdx(curr_idx)
            if curr_atom.GetAtomicNum() != 6:
                continue
            acyl_atoms.add(curr_idx)
            for neighbor in curr_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in acyl_atoms:
                    stack.append(neighbor.GetIdx())
        
        fatty_acyl_chain_length = len(acyl_atoms)
        # Check if the fatty acyl chain length (in terms of carbon count) falls within the expected range.
        if 13 <= fatty_acyl_chain_length <= 22:
            reason = (f"Found thioester linking a fatty acyl chain with {fatty_acyl_chain_length} carbons "
                      f"(including the carbonyl carbon) within the expected range (13-22) attached to a CoA moiety")
            return True, reason
        else:
            reasons.append(f"Fatty acyl chain comprised of {fatty_acyl_chain_length} carbons; expected between 13 and 22")
    
    reason_str = "; ".join(reasons) if reasons else "No valid long-chain fatty acyl-CoA structure found"
    return False, reason_str