"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3‑substituted propionyl‑CoA(4-)
Definition:
  An acyl‐CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups 
  of any 3‑substituted propionyl‑CoA; major species at pH 7.3.
Heuristic criteria used here:
  1. The molecule must be valid.
  2. Overall formal charge == -4.
  3. It must contain a thioester group (a [#6](=O)[S] fragment).
  4. It must contain a characteristic CoA moiety fragment (we search for a pantetheine core fragment):
       SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP
  5. It must have at least 4 deprotonated oxygens ([O-]).
  6. It must have exactly 3 carbonyl groups ([CX3](=O)).
  7. By “extracting” the acyl chain (the part attached to the thioester carbonyl)
     we require that it is a simple, linear hydrocarbon (only carbon atoms in a straight chain)
     with no branching.
If all conditions are met, the molecule is classified as a 3‑substituted propionyl‑CoA(4-).
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
    # In the match, index 0 is the carbonyl carbon and index 1 is the sulfur.
    thio_match = thioester_matches[0]
    carbonyl_idx = thio_match[0]
    sulfur_idx = thio_match[1]
    
    # 4. Check for the CoA moiety fragment.
    coa_smarts = "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety fragment not found"
    
    # 5. Check for at least 4 deprotonated oxygens ([O-]).
    deprot_oxy_count = sum(1 for atom in mol.GetAtoms()
                           if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if deprot_oxy_count < 4:
        return False, f"Found {deprot_oxy_count} deprotonated oxygen(s); expected at least 4 for CoA(4-)"
    
    # 6. Count all carbonyl groups ([CX3](=O)) in the molecule.
    carbonyl_smarts = "[CX3](=O)"
    carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) != 3:
        return False, f"Found {len(carbonyl_matches)} carbonyl group(s); expected exactly 3 (1 acyl + 2 in CoA fragment)"
    
    # 7. Extract and validate the acyl chain fragment.
    # The acyl chain should be the part attached to the thioester carbonyl (via a single C-C bond).
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    # Identify the carbon neighbor (the acyl chain start) that is not the carbonyl oxygen.
    acyl_neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors()
                      if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != sulfur_idx]
    if len(acyl_neighbors) != 1:
        return False, "Unexpected connectivity at thioester carbonyl; cannot isolate acyl chain"
    acyl_start = acyl_neighbors[0]
    
    # Instead of a full DFS that might wander into branches from the acyl chain,
    # we now attempt to “walk” a linear path along the acyl chain.
    chain_atom_indices = [acyl_start.GetIdx()]
    visited = set(chain_atom_indices)
    prev_idx = carbonyl_idx  # start from the carbonyl as 'previous'
    current_atom = acyl_start
    while True:
        # Get all neighboring carbon atoms, excluding the atom we came from.
        nbr_carbons = [nbr for nbr in current_atom.GetNeighbors()
                       if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_idx]
        if len(nbr_carbons) == 0:
            # Terminal carbon reached.
            break
        elif len(nbr_carbons) == 1:
            next_atom = nbr_carbons[0]
            # Check if we already visited this atom; if yes, break.
            if next_atom.GetIdx() in visited:
                break
            chain_atom_indices.append(next_atom.GetIdx())
            visited.add(next_atom.GetIdx())
            prev_idx = current_atom.GetIdx()
            current_atom = next_atom
        else:
            # Branch detected: more than one carbon neighbor.
            return False, "Acyl chain is branched; expected a simple fatty acyl chain with only carbon atoms"
    
    # Finally, verify that every atom in the acyl chain is carbon.
    for idx in chain_atom_indices:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6:
            return False, "Acyl chain contains non-carbon atoms; not a simple fatty acyl chain"
    
    # All conditions passed.
    return True, ("Molecule contains a thioester group, a CoA core fragment, the expected number of deprotonated oxygens, "
                  "exactly 3 carbonyl groups, overall charge -4 and an unbranched hydrocarbon-only acyl chain")

# Example usage (for testing purposes):
if __name__ == '__main__':
    # Test with one of the provided true positive SMILES strings.
    test_smiles = ("CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)"
                   "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)"
                   "OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    res, reason = is_3_substituted_propionyl_CoA_4__(test_smiles)
    print(res, reason)