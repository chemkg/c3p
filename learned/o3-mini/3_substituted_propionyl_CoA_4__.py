"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3‑substituted propionyl‑CoA(4-)
Definition:
  An acyl‐CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups 
  of any 3‑substituted propionyl‑CoA; major species at pH 7.3.
Heuristic criteria used here:
  1. The molecule must be a valid structure.
  2. It must contain a thioester group (a [#6](=O)[S] fragment).
  3. It must contain a characteristic CoA moiety fragment (we search for a pantetheine core):
       SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP
  4. It must have at least four deprotonated oxygens ([O-]).
  5. It must have exactly three carbonyl groups ([CX3](=O)).
  6. In addition, its overall formal charge must be –4.
  7. Finally, by “extracting” the acyl chain (the part attached to the thioester carbonyl),
     we require that this fragment contains only carbon atoms (i.e. it is a simple fatty acyl chain)
     as expected for a 3‑substituted propionyl‑CoA.
If all conditions are met the molecule is classified as a 3‑substituted propionyl‑CoA(4-).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule conforms to the 3‑substituted propionyl‑CoA(4-) criteria.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule matches the criteria, False otherwise.
        str: Explanation of the classification decision.
    """
    # 1. Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Check overall formal charge is exactly -4.
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -4:
        return False, f"Net formal charge is {net_charge}; expected -4"
    
    # 3. Check for a thioester group.
    # SMARTS: a carbonyl carbon ([#6]) double‐bonded to O and attached to S.
    thioester_smarts = "[#6](=O)[S]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (C(=O)S) group found"
    # Use the first match for further analysis.
    thio_match = thioester_matches[0]
    # In the match, index 0 is the carbonyl carbon and index 1 is the sulfur.
    carbonyl_idx = thio_match[0]
    sulfur_idx = thio_match[1]
    
    # 4. Check for the CoA moiety fragment.
    coa_smarts = "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety fragment not found"
    
    # 5. Check for at least four deprotonated oxygens ([O-]).
    deprot_oxy_count = sum(1 for atom in mol.GetAtoms()
                           if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if deprot_oxy_count < 4:
        return False, f"Found {deprot_oxy_count} deprotonated oxygen(s); expected at least 4 for CoA(4-)"
    
    # 6. Count all carbonyl groups [CX3](=O) in the molecule.
    carbonyl_smarts = "[CX3](=O)"
    carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    n_carbonyl = len(carbonyl_matches)
    if n_carbonyl != 3:
        return False, f"Found {n_carbonyl} carbonyl group(s); expected exactly 3 (1 acyl + 2 in CoA fragment)"
    
    # 7. Extract and validate the acyl chain fragment.
    # The acyl chain is defined as the part attached to the thioester carbonyl.
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    # Get neighbors of carbonyl carbon and pick the one that is not the carbonyl oxygen.
    acyl_neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() 
                      if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != sulfur_idx]
    if len(acyl_neighbors) != 1:
        return False, "Unexpected connectivity at thioester carbonyl; cannot isolate acyl chain"
    acyl_start = acyl_neighbors[0].GetIdx()
    
    # We perform a simple DFS from the acyl_start atom without crossing back to the carbonyl atom.
    # We also do not follow the bond to the sulfur (which leads to the CoA moiety).
    visited = set()
    stack = [acyl_start]
    acyl_atoms = set()
    exclude = {carbonyl_idx, sulfur_idx}
    while stack:
        aidx = stack.pop()
        if aidx in visited:
            continue
        visited.add(aidx)
        acyl_atoms.add(aidx)
        atom = mol.GetAtomWithIdx(aidx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in exclude or nbr_idx in visited:
                continue
            # Continue traversal only if the neighbor is connected via a single bond
            # (this is a heuristic to remain in the acyl chain region).
            if mol.GetBondBetweenAtoms(aidx, nbr_idx).GetBondType() in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
                stack.append(nbr_idx)
    
    # Now check that every atom in the acyl chain is carbon.
    for aidx in acyl_atoms:
        atom = mol.GetAtomWithIdx(aidx)
        if atom.GetAtomicNum() != 6:
            return False, "Acyl chain contains heteroatoms; not a simple fatty acyl chain as required for 3‑substituted propionyl‑CoA"
    
    # Passed all criteria.
    return True, ("Molecule contains a thioester group, a CoA core fragment, the expected number of deprotonated oxygens, "
                  "exactly 3 carbonyl groups, overall charge -4 and a hydrocarbon-only acyl chain")

# Example usage (for testing purposes):
if __name__ == '__main__':
    # Test with one of the provided true positive SMILES strings.
    test_smiles = ("CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)"
                   "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)"
                   "OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    res, reason = is_3_substituted_propionyl_CoA_4__(test_smiles)
    print(res, reason)