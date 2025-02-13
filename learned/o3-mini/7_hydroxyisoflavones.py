"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
#!/usr/bin/env python3
"""
Classifies: 7-hydroxyisoflavones
Definition: A hydroxyisoflavone compound having a hydroxy (–OH) group specifically at the 7–position of the isoflavone (3–phenylchromen-4-one) core.
Heuristic algorithm:
1. Look for a relaxed isoflavone core using a SMARTS pattern.
   (The SMARTS "c1ccc2c(c1)occ2=O" represents the fused bicyclic (chromen-4-one) system.)
2. Within the matched core, identify the carbonyl carbon (one that is double-bonded to an oxygen).
3. Iterate over the aromatic carbons in the core that bear an –OH substituent.
4. Evaluate whether one such carbon is exactly three bonds away (topologically) from the carbonyl carbon.
5. If so, classify as a 7-hydroxyisoflavone.
Note: This rule‐based procedure is heuristic and may not capture all isoflavone derivatives.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    The function requires that the molecule contains an isoflavone core and then
    checks for the presence of an aromatic hydroxyl group (–OH) on the core at a location
    that is three bonds away from the carbonyl group (heuristically corresponding to the 7–position).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 7-hydroxyisoflavone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a relaxed SMARTS for the isoflavone core (chromen-4-one)
    # This pattern represents a fused bicyclic aromatic system with a lactone (carbonyl) group.
    core_smarts = "c1ccc2c(c1)occ2=O"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in SMARTS pattern for isoflavone core."
    
    # Check if the molecule contains the core
    if not mol.HasSubstructMatch(core_query):
        return False, "Molecule does not contain the expected isoflavone core."
    
    core_matches = mol.GetSubstructMatches(core_query)
    # Use the first match found
    core_match = core_matches[0]
    core_match_set = set(core_match)
    
    # Within the core match, identify the carbonyl carbon.
    # We search for a carbon which is double-bonded to an oxygen.
    carbonyl_idx = None
    for idx in core_match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:  # carbon
            # Look at bonds to see if one is a double bond to oxygen
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() == 2:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8:
                        carbonyl_idx = idx
                        break
            if carbonyl_idx is not None:
                break
    if carbonyl_idx is None:
        return False, "Isoflavone core found but no carbonyl carbon was detected in the core."
    
    # Search within the core for an aromatic carbon that has an -OH group
    # and check if its shortest topological distance to the carbonyl carbon is 3 bonds.
    for idx in core_match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
            continue
        # Check if this aromatic carbon is substituted with a hydroxyl group.
        has_hydroxy = False
        for nbr in atom.GetNeighbors():
            # Look for oxygen atoms.
            if nbr.GetAtomicNum() != 8:
                continue
            # Check if the oxygen is a hydroxyl: typically, it has at least one hydrogen (even if implicit)
            # and is not involved in a double bond (to exclude carbonyl oxygen).
            if nbr.GetDegree() >= 1 and nbr.GetFormalCharge() == 0:
                # Get total H count (implicit + explicit)
                if nbr.GetTotalNumHs() > 0:
                    has_hydroxy = True
                    break
        if not has_hydroxy:
            continue
        
        # Compute the shortest topological distance from this aromatic carbon to the carbonyl carbon.
        path = rdmolops.GetShortestPath(mol, idx, carbonyl_idx)
        # The number of bonds is (number of atoms in path - 1)
        if len(path) - 1 == 3:
            return True, ("Molecule contains an isoflavone core and an aromatic -OH located "
                          "three bonds from the carbonyl (consistent with a 7-hydroxy substitution).")
    
    return False, ("Molecule contains an isoflavone core but no aromatic -OH was found at the expected "
                   "position (three bonds from the carbonyl).")

# Example (for local testing):
# if __name__ == "__main__":
#     test_smiles = "Oc1ccc2c(c1)occ(-c1ccccc1)c2=O"  # 7-hydroxyisoflavone
#     result, reason = is_7_hydroxyisoflavones(test_smiles)
#     print(result, reason)