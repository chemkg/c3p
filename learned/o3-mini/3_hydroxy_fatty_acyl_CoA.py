"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
Definition: A hydroxy fatty acyl-CoA that results from the formal condensation
            of the thiol group of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.

A valid molecule must contain:
 - A recognizable CoA moiety (using a typical pantetheine fragment pattern)
 - A thioester group (C(=O)S)
 - In the acyl chain linked at the thioester carbonyl, if the carbonyl is C1,
   the beta carbon (C3) must bear a neutral –OH group that is actually protonated.
   
Our strategy:
 1. Parse the SMILES and quit if invalid.
 2. Check for a CoA fragment (using a SMARTS) and record the indices that belong to it.
 3. Find a thioester group (using SMARTS "C(=O)S").
 4. For each thioester found, take the carbonyl carbon (C1). Then, among its neighbors,
    pick potential acyl-chain (alpha) carbons that are not part of the CoA fragment.
 5. For each alpha carbon (C2), look at its neighbors (excluding C1) for a beta carbon (C3).
 6. For each candidate beta carbon, check that it bears an oxygen by a single bond.
    In addition, require that the oxygen is “protonated” (total hydrogen count > 0)
    and has no formal charge.
 7. If any thioester-acyl side chain meets these criteria, classify as a 3-hydroxy fatty acyl-CoA.
    
If any step fails, return False with a reason.
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    
    The algorithm:
      1. Parses the SMILES.
      2. Checks for a Coenzyme A fragment (using a typical pantetheine-like SMARTS).
      3. Searches for a thioester group (“C(=O)S”).
      4. For each thioester, from the carbonyl carbon (C1), it finds a candidate α‐carbon (C2)
         that is not part of the CoA fragment.
      5. From the α‐carbon it finds a candidate β‐carbon (C3) and checks that C3 carries an -OH
         substituent by a single bond where the oxygen is neutral and has an attached hydrogen.
      6. If any such valid acyl chain is found, the molecule is classified as a 3‐hydroxy fatty acyl-CoA.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as 3-hydroxy fatty acyl-CoA, False otherwise.
      str: Reason for the classification decision.
    """
    # 1. Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Look for a Coenzyme A fragment.
    # Use a SMARTS pattern resembling the typical pantetheine unit in CoA.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_frag = Chem.MolFromSmarts(coa_smarts)
    coa_matches = mol.GetSubstructMatches(coa_frag)
    if not coa_matches:
        return False, "Coenzyme A fragment missing"
    # We combine all matching atom indices (if there are multiple CoA fragments, take the union)
    coa_atoms = set()
    for match in coa_matches:
        coa_atoms.update(match)
    
    # 3. Look for thioester groups via SMARTS: "C(=O)S"
    thioester_smarts = "C(=O)S"
    thioester_frag = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_frag)
    if not thioester_matches:
        return False, "No thioester group (acyl-CoA linkage) found"
    
    # 4. For each thioester group found, examine the acyl chain.
    # In our SMARTS match the first atom is the carbonyl carbon (C1),
    # the others are the C=O oxygen and the sulfur.
    valid_acyl = False
    for match in thioester_matches:
        if len(match) < 3:
            continue  # not a proper match, skip
        c1 = mol.GetAtomWithIdx(match[0])
        # 5. Identify alpha carbons: neighbors of the carbonyl carbon (C1)
        # Exclude atoms in the thioester match (O and S) and also exclude any atom that belongs to the CoA fragment.
        alpha_candidates = []
        for nbr in c1.GetNeighbors():
            if nbr.GetIdx() in match[1:]:
                continue
            if nbr.GetIdx() in coa_atoms:
                continue
            if nbr.GetAtomicNum() == 6:  # carbon
                alpha_candidates.append(nbr)
        if not alpha_candidates:
            continue  # no valid alpha found on this thioester
        
        # 6. For each alpha candidate (C2), look for beta carbons (C3)
        for alpha in alpha_candidates:
            beta_candidates = []
            for nbr in alpha.GetNeighbors():
                if nbr.GetIdx() == c1.GetIdx():  # avoid going back to carbonyl
                    continue
                if nbr.GetAtomicNum() == 6:  # must be a carbon
                    beta_candidates.append(nbr)
            if not beta_candidates:
                continue
            # For each beta candidate, check if it carries an -OH substituent
            for beta in beta_candidates:
                # Look over beta's neighbors
                for beta_nbr in beta.GetNeighbors():
                    # Identify an oxygen attached by a single bond.
                    if beta_nbr.GetAtomicNum() != 8:
                        continue
                    bond = mol.GetBondBetweenAtoms(beta.GetIdx(), beta_nbr.GetIdx())
                    if bond is None or bond.GetBondTypeAsDouble() != 1:
                        continue
                    # Require that the oxygen is not deprotonated
                    if beta_nbr.GetFormalCharge() != 0:
                        continue
                    # Check that the oxygen appears to have at least one hydrogen.
                    if beta_nbr.GetTotalNumHs() > 0:
                        valid_acyl = True
                        break
                if valid_acyl:
                    break
            if valid_acyl:
                break
        if valid_acyl:
            break
    
    if not valid_acyl:
        return False, "Thioester found but no acyl chain with a neutral -OH on the beta carbon (C3) detected"
    
    return True, "Molecule contains a CoA moiety, a thioester group, and an acyl chain with a protonated -OH at the beta position"

# (Optional) If run as a script, try some examples.
if __name__ == '__main__':
    # A few examples from the provided true positive list:
    test_smiles = [
        "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
        "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C1CCCCC1O",
    ]
    for s in test_smiles:
        res, reason = is_3_hydroxy_fatty_acyl_CoA(s)
        print(res, reason)