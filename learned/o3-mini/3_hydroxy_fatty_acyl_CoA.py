"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
Definition: A hydroxy fatty acyl-CoA that results from the formal condensation
            of the thiol group of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.
A valid molecule must contain:
 - A recognizable CoA moiety (using a SMARTS for the typical pantetheine fragment)
 - A thioester group (C(=O)S)
 - For at least one thioester linkage, the acyl chain (attached at the carbonyl)
   should be “fatty-acid–like”: that is, it must be a relatively unbranched chain of several carbons
   where, if the carbonyl is C1 and its alpha carbon is C2, then the beta carbon (C3) bears a neutral –OH.
   
Our revised strategy:
 1. Parse the SMILES.
 2. Check for a Coenzyme A substructure.
 3. From all thioester fragments (SMARTS "C(=O)S"), for each:
      a. Let c1 = carbonyl carbon.
      b. Identify alpha candidate (C2): a carbon neighbor of c1 that is not one of the matched thioester atoms and not part of the CoA fragment.
      c. Check that C2 is “linear” (i.e. has exactly one carbon neighbor besides c1). That neighbor becomes our beta candidate (C3).
      d. Verify that C3 carries an -OH substituent (an oxygen attached via a single bond that is neutral and shows at least one hydrogen).
      e. In addition, perform a DFS starting at C2 (excluding C1) to estimate the acyl chain length. We require a minimum length (here ≥6 carbons) to rule out small or branched fragments.
 4. If any thioester linkage meets the above criteria, classify as a 3-hydroxy fatty acyl-CoA.
 Otherwise, return False with an appropriate reason.
"""

from rdkit import Chem

def get_max_chain_length(atom, parent_idx, visited):
    """
    Recursively find the maximum number of carbon atoms (counting the current)
    in a linear acyl chain starting from 'atom' (avoiding going back to 'parent_idx').
    'visited' is used to prevent cycles.
    """
    # Only consider carbon atoms
    if atom.GetAtomicNum() != 6:
        return 0
    max_length = 1  # count current carbon
    visited.add(atom.GetIdx())
    for nbr in atom.GetNeighbors():
        if nbr.GetIdx() == parent_idx:
            continue
        if nbr.GetIdx() in visited:
            continue
        if nbr.GetAtomicNum() == 6:
            # Continue down this bond
            length = 1 + get_max_chain_length(nbr, atom.GetIdx(), visited.copy())
            if length > max_length:
                max_length = length
    return max_length

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.

    The algorithm requires that:
      - The molecule contains a Coenzyme A moiety (identified by a pantetheine-like pattern).
      - It contains a thioester group (C(=O)S).
      - For at least one thioester, the acyl chain attached (starting at the carbonyl carbon)
        is linear and long enough (≥6 carbons in the chain) and its beta carbon (the second carbon
        of the acyl chain, counting the carbonyl as C1 and the alpha as C2) bears an -OH that is protonated.
        
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule meets the criteria, False otherwise.
      str: Explanation for the classification decision.
    """
    # 1. Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Check for a CoA fragment.
    # We use a SMARTS pattern that covers part of the pantetheine unit.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_frag = Chem.MolFromSmarts(coa_smarts)
    coa_matches = mol.GetSubstructMatches(coa_frag)
    if not coa_matches:
        return False, "Coenzyme A fragment missing"
    # Combine all matching atom indices from CoA fragment(s)
    coa_atoms = set()
    for match in coa_matches:
        coa_atoms.update(match)
    
    # 3. Look for thioester groups (SMARTS "C(=O)S")
    thioester_smarts = "C(=O)S"
    thioester_frag = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_frag)
    if not thioester_matches:
        return False, "No thioester group (acyl-CoA linkage) found"
    
    # 4. For each thioester match, try to verify the acyl chain connectivity and the beta -OH.
    for match in thioester_matches:
        # In the "C(=O)S" SMARTS match:
        # match[0] is the carbonyl carbon (C1),
        # match[1] is the carbonyl oxygen, and
        # match[2] is the sulfur.
        if len(match) < 3:
            continue  # malformed match
        
        c1 = mol.GetAtomWithIdx(match[0])
        
        # 4a. Identify an alpha candidate (C2): a neighbor of c1 that is
        # a carbon (atomic number 6), is not the oxygen or sulfur from the match,
        # and is not part of the CoA fragment.
        alpha_candidates = []
        for nbr in c1.GetNeighbors():
            if nbr.GetIdx() in match[1:]:
                continue  # skip the O and S atoms in the thioester
            if nbr.GetIdx() in coa_atoms:
                continue  # skip atoms from the CoA fragment
            if nbr.GetAtomicNum() == 6:
                alpha_candidates.append(nbr)
        if not alpha_candidates:
            continue  # try next thioester if no eligible alpha carbon found
        
        for alpha in alpha_candidates:
            # 4b. For a typical fatty acid chain the alpha carbon (C2) should be linear --
            # that is, apart from its bond back to c1, it should have exactly one carbon neighbor.
            carbon_nbrs = [n for n in alpha.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != c1.GetIdx()]
            if len(carbon_nbrs) != 1:
                continue  # skip if alpha is branched
            
            beta = carbon_nbrs[0]  # candidate beta carbon (C3)
            
            # 4c. Check that beta carries an -OH group (oxygen attached via a single bond,
            # with formal charge 0 and at least one hydrogen).
            oh_found = False
            for beta_nbr in beta.GetNeighbors():
                if beta_nbr.GetAtomicNum() != 8:
                    continue
                bond = mol.GetBondBetweenAtoms(beta.GetIdx(), beta_nbr.GetIdx())
                if bond is None or bond.GetBondTypeAsDouble() != 1:
                    continue
                if beta_nbr.GetFormalCharge() != 0:
                    continue
                if beta_nbr.GetTotalNumHs() < 1:
                    continue
                oh_found = True
                break
            if not oh_found:
                continue  # no proper beta -OH on this chain; try next alpha
            
            # 4d. Check the acyl chain length using a DFS starting at alpha.
            # We count the maximum number of consecutive carbon atoms (excluding c1).
            chain_length = get_max_chain_length(alpha, c1.GetIdx(), set())
            if chain_length < 6:
                # Too short to be considered a fatty acyl chain
                continue
            
            # If all criteria are met, then we classify this molecule as a valid 3-hydroxy fatty acyl-CoA.
            return True, ("Molecule contains a CoA moiety, a thioester group and an acyl chain "
                          "with adequate length having a protonated -OH at the beta (C3) position.")
    
    # If we never found a valid thioester-acyl chain meeting criteria:
    return False, ("Thioester linkage found but no acyl chain meeting criteria: "
                   "either lacking a linear fatty-acid chain (≥6 carbons) or missing a protonated -OH on the beta carbon.")


# (Optional) Main block to try a few examples.
if __name__ == '__main__':
    # Some examples from the provided outcomes:
    test_smiles = [
        "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
        "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C1CCCCC1O",  # a true positive example
        # A false positive example (e.g. 3-hydroxyisovaleryl-CoA shown with (4-) in the name)
        "CC(C)(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12",
    ]
    for s in test_smiles:
        res, reason = is_3_hydroxy_fatty_acyl_CoA(s)
        print(res, reason)