"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: CHEBI:3-oxo-Δ(1) steroid
Definition: Any 3-oxo steroid that contains a double bond between positions 1 and 2.
Heuristic improvements:
  - Instead of using a single SMARTS pattern, iterate over candidate carbonyl carbons (ketone groups)
    in six-membered rings.
  - Check that the candidate carbonyl carbon (presumably the 3-oxo group) is indeed a ketone 
    (i.e. bonded via a double bond to oxygen and to two carbon atoms).
  - Then, for each candidate, look for an adjacent ring carbon (in the same six‐membered ring)
    that is part of a non-carbonyl double bond (a candidate for the Δ(1) double bond, usually between C1 and C2).
  - Also require that the molecule has at least 3 rings (typically steroids have 4 fused rings)
    and that the candidate carbon is part of a fused ring system.
Note: This heuristic does not cover all edge cases.
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Δ(1) steroid based on its SMILES string.
    
    Heuristic criteria:
      1. The molecule must have at least 3 rings (steroid nucleus is typically 4 fused rings).
      2. In at least one six-membered ring, there is a ketone carbon (a carbon with a C=O double bond,
         and with at least 2 carbon neighbours, so that it is a ketone not an aldehyde).
      3. One of the ring-neighbours of that ketone carbon participates in a non-carbonyl C=C double bond.
         (This is used as a proxy for the Δ(1) double bond between positions 1 and 2.)
      4. The candidate ketone should be in a fused ring system (i.e. be part of >1 ring).
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if the molecule qualifies as a 3-oxo-Δ(1) steroid, else False.
       str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    # Require at least 3 rings (steroids usually have 4 fused rings)
    if ring_info.NumRings() < 3:
        return False, "Molecule does not have enough rings (<3) to be a steroid"
    
    # Iterate over atoms to find candidate ketone carbons in six-membered rings.
    for atom in mol.GetAtoms():
        # Look for carbon atoms in a six-membered ring that might be a ketone carbon.
        if atom.GetSymbol() != "C":
            continue
        # RDKit provides a helper: IsInRingSize(n)
        try:
            if not atom.IsInRingSize(6):
                continue
        except AttributeError:
            # if IsInRingSize is not available, alternate approach using ring_info.NumAtomRings(idx)
            if ring_info.NumAtomRings(atom.GetIdx()) == 0:
                continue
        
        # Check for carbonyl: an sp2 carbon bonded via a double bond to an oxygen.
        carbonyl_found = False
        o_neighbor = None
        for bond in atom.GetBonds():
            # Check if the bond is double and the other atom is oxygen.
            if bond.GetBondTypeAsDouble() == 2:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetSymbol() == "O":
                    carbonyl_found = True
                    o_neighbor = nbr
                    break
        if not carbonyl_found:
            continue
        
        # Verify that the atom is a ketone (attached to at least two carbon atoms besides the oxygen)
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) < 2:
            continue  # might be an aldehyde
        
        # Check if this candidate is part of a fused ring system:
        # At least one atom in the candidate's rings should be in more than one ring.
        # We can check the candidate atom itself:
        if ring_info.NumAtomRings(atom.GetIdx()) < 2:
            # Optionally, one could check the ring partner atoms; here we require the ketone center
            # to be fused.
            continue
        
        # Now, look among the ring neighbors for a candidate double bond (non-carbonyl double bond)
        found_delta = False
        for nbr in atom.GetNeighbors():
            # We focus on neighbors that are carbons and are in a six-membered ring.
            if nbr.GetSymbol() != "C":
                continue
            try:
                if not nbr.IsInRingSize(6):
                    continue
            except AttributeError:
                continue
            # For the neighbor, check if any of its bonds (except the one linking back to our ketone)
            # is a double bond between two carbons (and not involved with oxygen).
            for bond in nbr.GetBonds():
                # Skip the bond that connects to the candidate ketone carbon.
                if bond.GetOtherAtom(nbr).GetIdx() == atom.GetIdx():
                    continue
                # Look for a double bond between carbons.
                if bond.GetBondTypeAsDouble() == 2:
                    other = bond.GetOtherAtom(nbr)
                    if other.GetSymbol() == "C":
                        # Additionally, ensure that this double bond is within a ring.
                        if bond.IsInRing():
                            found_delta = True
                            break
            if found_delta:
                break
        
        if found_delta:
            return True, ("Molecule contains a six-membered ring ketone (likely the 3-oxo group) with an adjacent "
                          "non-carbonyl C=C double bond (proxy for the Δ(1) bond) in a fused ring system.")
    
    return False, "No candidate six-membered ring ketone with adjacent non-carbonyl double bond (Δ(1) candidate) found"


# For testing purposes (you can remove or comment out the test cases)
if __name__ == "__main__":
    test_smiles = [
        # paraminabeolide B (expected True)
        "C[C@@H]([C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@H]3CC[C@]12COC(C)=O",
        # helvolic acid methyl ester (expected True)
        "C=1[C@@]2([C@@]3(CC[C@@]/4([C@@]([C@]3(C([C@H]([C@]2([C@@H](C(C1)=O)C)[H])OC(=O)C)=O)C)(C[C@@H](\\C4=C(\\CCC=C(C)C)/C(=O)OC)OC(=O)C)C)[H])[H])C",
        # estra-1,5(10)-diene-3,4,17-trione (expected True)
        "C[C@]12CC[C@H]3[C@@H](CCC4=C3C=CC(=O)C4=O)[C@@H]1CCC2=O",
    ]
    for smi in test_smiles:
        res, reason = is_3_oxo_Delta_1__steroid(smi)
        print("SMILES:", smi)
        print("Classification:", res)
        print("Explanation:", reason)
        print("----------")