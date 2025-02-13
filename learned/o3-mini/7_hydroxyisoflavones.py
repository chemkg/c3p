"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
#!/usr/bin/env python3
"""
Classifies: 7-hydroxyisoflavones
Definition: A hydroxyisoflavone compound having a hydroxy (–OH) group specifically at the 7–position of the isoflavone core.
Heuristic algorithm:
1. Verify that the molecule contains an isoflavone core (a 3–phenylchromen–4–one) using a SMARTS pattern.
2. Locate the carbonyl carbon (of the chromen-4-one portion) within that core.
3. Look for fused aromatic rings that are parts of the isoflavone core (the benzene rings in the fused bicyclic system).
4. Among those rings, search for an aromatic carbon bearing an –OH group whose topological distance to the carbonyl is 3.
   (In many isoflavones the 7–OH lies three bonds away from the carbonyl carbon.)
5. If such an –OH is found, classify the molecule as a 7–hydroxyisoflavone.
Note: This rule‐based procedure is heuristic and may miss edge cases or mis‐classify some molecules.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    The function requires that the molecule contains an isoflavone core, and then it
    checks for a hydroxyl group on the fused aromatic A ring that is three bonds away 
    (topologically) from the carbonyl carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 7-hydroxyisoflavone, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Step 1. Check for isoflavone core.
    # This SMARTS represents a 3-phenylchromen-4-one (isoflavone) core.
    core_smarts = "c1ccc2c(c1)oc(=O)c3ccccc23"
    core = Chem.MolFromSmarts(core_smarts)
    if core is None:
        return False, "Error in SMARTS pattern for isoflavone core."
    
    if not mol.HasSubstructMatch(core):
        return False, "Molecule does not contain the expected isoflavone core."
    
    # Get one match for the core.
    core_matches = mol.GetSubstructMatches(core)
    core_match = core_matches[0]  # use first match
    core_match_set = set(core_match)
    
    # Step 2. Find the carbonyl carbon in the core.
    # We look for a pattern matching a carbon with a double-bonded oxygen: [C]=O
    carbonyl_smarts = "[C](=O)"
    carbonyl_query = Chem.MolFromSmarts(carbonyl_smarts)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_query)
    carbonyl_core = None
    for match in carbonyl_matches:
        # If the carbonyl carbon is part of the core match, take it.
        if match[0] in core_match_set:
            carbonyl_core = match[0]
            break
    if carbonyl_core is None:
        return False, "Isoflavone core found but no carbonyl carbon detected."

    # Step 3. Get the ring information and count ring memberships.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    
    # Step 4. Look for aromatic six-membered rings that are part of the core match.
    candidate_rings = []
    for ring in rings:
        # Consider 6-membered aromatic rings fully contained in the isoflavone core.
        ring_set = set(ring)
        if len(ring) != 6:
            continue
        if not ring_set.issubset(core_match_set):
            continue
        # Check that all atoms in this ring are aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        candidate_rings.append(ring)
    
    if not candidate_rings:
        return False, "Isoflavone core found but no fused aromatic ring (A ring) identified within the core."
    
    # Step 5. Now, for each candidate ring, search for an aromatic carbon carrying -OH
    # and check if its shortest topological path to the carbonyl carbon is 3 bonds.
    for ring in candidate_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look only at aromatic carbons.
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                continue
            # Check neighbors for -OH (oxygen with explicit hydrogen(s))
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    # Found a hydroxyl on this aromatic carbon.
                    # Now, compute the shortest path length between this carbon and the carbonyl carbon.
                    path = rdmolops.GetShortestPath(mol, idx, carbonyl_core)
                    if len(path) - 1 == 3:
                        return True, "Molecule contains an isoflavone core and an aromatic -OH three bonds away from the carbonyl (consistent with a 7-OH)."
    
    return False, "Molecule contains an isoflavone core but no hydroxyl group at the expected (7-position) location was found."


# Example test (uncomment for local testing):
# if __name__ == "__main__":
#     test_smiles = "Oc1ccc2c(c1)occ(-c1ccccc1)c2=O"  # 7-hydroxyisoflavone
#     result, reason = is_7_hydroxyisoflavones(test_smiles)
#     print(result, reason)