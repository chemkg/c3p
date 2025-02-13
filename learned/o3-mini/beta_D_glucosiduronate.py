"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate
A carbohydrate acid derivative anion obtained by deprotonation of the carboxy group 
of any beta-D-glucuronic acid; major species at pH 7.3.
In glucuronic acid the primary alcohol (CH2OH) of glucose is oxidized to a carboxylate 
group (C(=O)[O-]) and the aglycone is attached via the oxygen on the anomeric carbon.
This improved heuristic first finds a carboxylate group and then verifies that the rest of 
the moiety is part of a six-membered (pyranose) ring (with one ring oxygen).
"""

from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule contains a beta-D-glucosiduronate moiety based on its SMILES string.
    
    The method uses a two-step heuristic:
      1. It searches for a carboxylate group, defined as C(=O)[O-].
      2. For each found carboxylate group, it examines the carbon atom attached
         to the carboxylate (other than the oxygens of the carboxylate). This neighbor,
         if present, should be part of a six-membered ring (a pyranose) in which there is exactly 
         one ring oxygen. Such a pattern is typical of glucuronic acid and its derivatives.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if a beta-D-glucosiduronate moiety is found, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for a deprotonated carboxylate group.
    carboxylate_smarts = "C(=O)[O-]"
    carboxylate_query = Chem.MolFromSmarts(carboxylate_smarts)
    if carboxylate_query is None:
        return False, "Error in SMARTS pattern definition"
    
    matches = mol.GetSubstructMatches(carboxylate_query)
    if not matches:
        return False, "No carboxylate group found in the molecule"
    
    # Get ring information once.
    ring_info = mol.GetRingInfo()
    
    # Iterate over each carboxylate match.
    for match in matches:
        # In the match tuple, match[0] is the carboxylate carbon and match[1] and match[2]
        # are the oxygens (one double-bonded and one carrying a negative charge).
        carboxylate_c = mol.GetAtomWithIdx(match[0])
        # Identify the neighbor that is not part of the C=O group.
        neighbor_indices = [nbr.GetIdx() for nbr in carboxylate_c.GetNeighbors()]
        # Exclude the oxygen atoms that make up the carboxylate.
        candidate_neighbors = [idx for idx in neighbor_indices if idx not in (match[1], match[2])]
        if not candidate_neighbors:
            # No neighbor found (should not normally occur).
            continue
        # Usually there is one neighbor (the carbon to which the carboxylate is attached).
        for nbr_idx in candidate_neighbors:
            nbr_atom = mol.GetAtomWithIdx(nbr_idx)
            if not nbr_atom.IsInRing():
                continue  # This candidate is not in any ring.
            # Examine all rings that include this atom.
            for ring in ring_info.AtomRings():
                if nbr_idx in ring and len(ring) == 6:
                    # Count how many atoms in this ring are oxygen.
                    ring_oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                    if ring_oxygen_count == 1:
                        # We found a carboxylate attached to a six-membered ring with one oxygen.
                        return True, ("Contains beta-D-glucosiduronate moiety: a carboxylate (C(=O)[O-]) "
                                      "attached to a pyranose ring with one ring oxygen.")
    # If no candidate ring environment is found:
    return False, "Beta-D-glucosiduronate moiety not found in the molecule"

# Example usage (uncomment to test):
# smiles_example = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3([C@@H](C[C@@]2(C[C@@H](C1)O)[H])O)[H])(CC[C@@]4([C@@H](CCC(O[C@@H]5O[C@@H]([C@H]([C@@H]([C@H]5O)O)O)C([O-])=O)=O)C)[H])[H])C)[H])C"
# result, reason = is_beta_D_glucosiduronate(smiles_example)
# print(result, reason)