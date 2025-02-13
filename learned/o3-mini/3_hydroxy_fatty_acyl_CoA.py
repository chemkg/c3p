"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
Definition: A hydroxy fatty acyl-CoA that results from the formal condensation 
            of the thiol group of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.

A molecule in this class must contain:
 - A recognizable CoA moiety (the pantetheine fragment, e.g. SCCNC(=O)CCNC(=O) is typical)
 - A thioester group (C(=O)S) linking an acyl chain to CoA.
 - In the acyl chain attached via the thioester, if one numbers the carbonyl carbon as C1,
   then the beta carbon (C3) must bear a hydroxyl (-OH) substituent.
   
Our strategy:
 1. Parse the SMILES.
 2. Look for a CoA fragment using a SMARTS query.
 3. Look for thioester groups (C(=O)S).
 4. For each thioester, from the carbonyl carbon (C1) find all carbon neighbors (the α‐carbons, C2),
    then from each α‐carbon find further carbon neighbors (the β‐carbons, C3). For each β,
    check if it is substituted by –OH (by looking for an oxygen attached by a single bond that has at least one hydrogen).
 5. If at least one valid acyl chain is found, classify as 3-hydroxy fatty acyl-CoA.
 
If parsing fails or if no match is found, the function returns False with a reason.
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    
    The algorithm:
      1. Parses the SMILES string.
      2. Checks for a Coenzyme A fragment by matching a characteristic SMARTS.
      3. Searches for thioester groups (C(=O)S) in the molecule.
      4. For each thioester found, it:
           a. Identifies the carbonyl carbon (C1).
           b. Finds candidate α‐carbons (neighbors of C1, excluding the oxygen in C=O and the sulfur).
           c. For each α‐carbon, finds candidate β‐carbons (neighbors of α, excluding the path back to C1).
           d. Checks if any β‐carbon carries an –OH substituent: that is, an oxygen attached via a single bond
              and the oxygen atom has an implicit hydrogen (or total H > 0).
      5. If such a pattern is found, the molecule is classified as a 3-hydroxy fatty acyl-CoA.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as 3-hydroxy fatty acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check presence of a Coenzyme A fragment.
    # Here we use a fragment typical for many CoA derivatives – adjust if necessary.
    cof_smarts = "SCCNC(=O)CCNC(=O)"
    coa_frag = Chem.MolFromSmarts(cof_smarts)
    if not mol.HasSubstructMatch(coa_frag):
        return False, "Coenzyme A fragment missing"
    
    # Look for thioester groups with SMARTS "C(=O)S"
    thioester_smarts = "C(=O)S"
    thioester_frag = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_frag)
    if not thioester_matches:
        return False, "No thioester group (acyl-CoA linkage) found"
    
    # For each thioester group found, examine the acyl chain
    # Each match returns a tuple: (carbonyl carbon, carbonyl oxygen, sulfur)
    valid_acyl_found = False
    for match in thioester_matches:
        # Define the carbonyl carbon (C1) as the first atom in the match.
        c1 = mol.GetAtomWithIdx(match[0])
        
        # Get neighbors of C1 that are carbons and are not the oxygen (double-bonded) or sulfur.
        alpha_candidates = []
        for nbr in c1.GetNeighbors():
            # Exclude if neighbor index is one of the oxygen or sulfur from the thioester pattern.
            if nbr.GetIdx() in match[1:]:
                continue
            if nbr.GetAtomicNum() == 6:  # Carbon
                alpha_candidates.append(nbr)
        if not alpha_candidates:
            continue  # No possible alpha carbon on this thioester
        
        # For each candidate alpha carbon (C2) try to find a beta carbon (C3) with an -OH.
        for alpha in alpha_candidates:
            # From alpha, get neighbors that are carbons and are not back to C1.
            beta_candidates = []
            for nbr in alpha.GetNeighbors():
                if nbr.GetIdx() == c1.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:  # keep carbons only
                    beta_candidates.append(nbr)
            if not beta_candidates:
                continue
            
            # For each beta candidate, search for an –OH substituent.
            for beta in beta_candidates:
                # Look at all neighbors of the beta carbon.
                for nbr in beta.GetNeighbors():
                    # Check if the neighbor is oxygen.
                    if nbr.GetAtomicNum() != 8:
                        continue
                    # Check that the bond is a single bond.
                    bond = mol.GetBondBetweenAtoms(beta.GetIdx(), nbr.GetIdx())
                    if bond is None or bond.GetBondTypeAsDouble() != 1:
                        continue
                    # Check if this oxygen has at least one hydrogen (the hydrogen count may be implicit).
                    if nbr.GetTotalNumHs() > 0:
                        valid_acyl_found = True
                        break
                if valid_acyl_found:
                    break
            if valid_acyl_found:
                break
        if valid_acyl_found:
            break

    if not valid_acyl_found:
        return False, "Thioester found but no acyl chain with a hydroxyl on the beta carbon (C3)"
    
    return True, "Molecule contains a CoA moiety, a thioester group, and an acyl chain with an -OH on the beta carbon"

# (Optional) For testing:
if __name__ == '__main__':
    # Example SMILES from the provided true positive list:
    smiles_examples = [
        "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
        "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C1CCCCC1O"  # 2-hydroxycyclohexane-1-carbonyl-CoA (false negative previously)
    ]
    for s in smiles_examples:
        result, reason = is_3_hydroxy_fatty_acyl_CoA(s)
        print(result, reason)