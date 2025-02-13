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
      1. A CoA moiety defined by the presence of a pantetheine fragment (SCCNC(=O)CCNC(=O)).
      2. A thioester group ([CX3](=O)[SX2]) linking the fatty acyl chain to the CoA.
      3. A fatty acyl chain (including the carbonyl carbon) that contains between 13 and 22 carbon atoms.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a long-chain fatty acyl-CoA, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Detect CoA moiety.
    # Instead of looking for an adenine fragment, we search for a pantetheine fragment characteristic for CoA:
    # The fragment "SCCNC(=O)CCNC(=O)" appears in many acyl-CoA molecules.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_frag = Chem.MolFromSmarts(coa_smarts)
    if coa_frag is None:
        return False, "Invalid CoA SMARTS pattern"
    if not mol.HasSubstructMatch(coa_frag):
        return False, "No CoA (coenzyme A) moiety detected (missing pantetheine fragment)"
    
    # 2. Detect the thioester group linking the fatty acid chain to CoA.
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if thioester is None:
        return False, "Invalid thioester SMARTS pattern"
    
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester group (linking fatty acyl chain to CoA) detected"
    
    reasons = []
    valid_match_found = False
    
    # 3. For each thioester group, try to isolate the fatty acyl chain.
    for match in thioester_matches:
        # match returns indices for the carbonyl carbon and the sulfur.
        carbonyl_idx, sulfur_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the neighbor of the carbonyl that is not the sulfur;
        # this is assumed to be the start of the fatty acyl chain.
        acyl_neighbors = []
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetIdx() == sulfur_idx:
                continue
            # Exclude oxygen (likely the double-bonded oxygen)
            if neighbor.GetAtomicNum() == 6:
                acyl_neighbors.append(neighbor)
        if len(acyl_neighbors) != 1:
            reasons.append("Ambiguous or missing acyl chain for one thioester match")
            continue
        acyl_start_idx = acyl_neighbors[0].GetIdx()
        
        # Now perform a DFS starting at the acyl chain neighbor.
        # We count carbon atoms reachable via carbon-carbon bonds (ignoring heteroatoms)
        # and we include the carbonyl carbon in the count.
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
            # Traverse only carbon-carbon bonds.
            for nb in curr_atom.GetNeighbors():
                if nb.GetAtomicNum() == 6 and nb.GetIdx() not in acyl_atoms:
                    stack.append(nb.GetIdx())
        
        fatty_acyl_chain_length = len(acyl_atoms)
        if 13 <= fatty_acyl_chain_length <= 22:
            valid_match_found = True
            reason = (f"Found thioester linking a fatty acyl chain with {fatty_acyl_chain_length} carbons "
                      f"(including the carbonyl carbon) within expected range (13-22) attached to a CoA moiety")
            return True, reason
        else:
            reasons.append(f"Fatty acyl chain comprised of {fatty_acyl_chain_length} carbons; expected between 13 and 22")
    
    if not valid_match_found:
        reason_str = "; ".join(reasons) if reasons else "No appropriate thioester fatty acyl chain found"
        return False, reason_str

    return False, "No valid long-chain fatty acyl-CoA structure found"