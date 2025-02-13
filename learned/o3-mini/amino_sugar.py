"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies amino sugars.
Definition: An amino sugar is a sugar molecule in which one or more alcoholic hydroxyl groups 
             attached to the sugar‐ring have been replaced by an unsubstituted amino group.
             
Improved algorithm:
  1. Parse the SMILES string.
  2. Identify candidate rings by:
       • Considering only rings of size 5 or 6.
       • Requiring that exactly one atom in the ring is oxygen and the remaining are carbons.
       • Rejecting rings that share atoms with any other ring (i.e. fused rings).
  3. For each candidate ring, look at exocyclic substituents off the ring carbons.
       • Only consider substituents that are a one‐atom branch (i.e. not more than a single atom)
         that are not part of any ring.
       • Count as “hydroxyl” if the neighbor is oxygen and has at least one explicit hydrogen.
       • Count as “free amino” only if the neighbor is nitrogen with exactly two explicit hydrogens
         and with no evidence of a C=O bond that might indicate an amide.
       • Require that there are at least 3 polar substituents on the sugar ring and at least one of them is free –NH2.
  4. If at least one candidate ring meets these criteria, return True.
  
If no candidate ring is found, or no candidate has a free amino substituent, return False.
"""

from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an amino sugar, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We update ring perception.
    mol.UpdatePropertyCache(strict=False)
    ring_info = mol.GetRingInfo()
    candidate_rings = []
    
    # Check each ring reported
    for ring in ring_info.AtomRings():
        # Accept only size 5 or 6 rings
        if len(ring) not in [5, 6]:
            continue
        # Get atoms in the ring
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check that all atoms are not aromatic (we want non‐aromatic sugar rings)
        if any(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        # Require exactly one oxygen and rest carbon
        oxygen_count = 0
        all_carbon = True
        for atom in ring_atoms:
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            elif atom.GetAtomicNum() != 6:
                all_carbon = False
                break
        if not all_carbon or oxygen_count != 1:
            continue
        # Check that none of the ring atoms are in a second ring (i.e. isolated candidate)
        fused = False
        for idx in ring:
            if len(ring_info.AtomRings(idx)) > 1:
                fused = True
                break
        if fused:
            continue
        # Candidate ring passes so far – add as a set of indices for simple membership testing
        candidate_rings.append(set(ring))
    
    if not candidate_rings:
        return False, "No isolated sugar-like ring (5- or 6-membered, non-aromatic with exactly one oxygen) detected."
    
    # Helper: check if a neighbor atom qualifies as a hydroxyl (–OH) substituent.
    def is_hydroxyl(neighbor):
        if neighbor.GetAtomicNum() != 8:
            return False
        if neighbor.IsInRing():
            return False
        # Count explicit hydrogens (using GetTotalNumHs with explicit only)
        # (This simple check assumes at least one hydrogen on the oxygen to form -OH)
        if neighbor.GetTotalNumHs(includeNeighbors=False) < 1:
            return False
        return True

    # Helper: check if a neighbor atom qualifies as a free amino (–NH2) substituent.
    def is_free_amino(neighbor):
        if neighbor.GetAtomicNum() != 7:
            return False
        if neighbor.IsInRing():
            return False
        # Require exactly 2 explicit hydrogens.
        if neighbor.GetTotalNumHs(includeNeighbors=False) != 2:
            return False
        # Check bonds from the nitrogen: if any double bond to oxygen (which may indicate an amide) then fail.
        for bond in neighbor.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0:
                other = bond.GetOtherAtom(neighbor)
                if other.GetAtomicNum() == 8:
                    return False
        return True

    # Now, for each candidate ring, look at substituents attached to ring carbons.
    for ring_set in candidate_rings:
        polar_substituent_count = 0
        found_amino = False
        # Only consider ring atoms that are carbon (ignoring the ring oxygen)
        for idx in ring_set:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Loop over neighbors that are not in the same ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                # Only consider substituents that are single-atom (to avoid long chains)
                # (One way is to demand that the neighbor has degree 1 or, if degree>1, then all bonds except the one
                # from the ring lead only to H's.)
                if nbr.GetDegree() > 1:
                    # Count non-hydrogen heavy-atom neighbors (other than the ring atom)
                    heavy = 0
                    for nn in nbr.GetNeighbors():
                        if nn.GetIdx() == atom.GetIdx():
                            continue
                        if nn.GetAtomicNum() != 1:
                            heavy += 1
                    if heavy > 0:
                        continue
                # Test for hydroxyl or free amino substituents
                if is_hydroxyl(nbr):
                    polar_substituent_count += 1
                elif is_free_amino(nbr):
                    polar_substituent_count += 1
                    found_amino = True
        # Require that the candidate sugar ring bears at least three small polar substituents 
        # AND that at least one is a free amino group.
        if found_amino and polar_substituent_count >= 3:
            return True, ("Molecule contains an isolated sugar-like ring (5- or 6-membered with one oxygen) " +
                          "decorated with at least three small polar substituents (including a free -NH2) " +
                          "consistent with an amino sugar structure.")
    
    return False, "No candidate sugar ring with a free -NH2 replacing an -OH group was detected."

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "N[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O",  # 3-amino-3-deoxy-D-glucopyranose, expected True
        "O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H](O)[C@H]1O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)CO[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O[C@@H]7O[C@@H]([C@@H](O[C@@H]8O[C@@H]([C@H](O)[C@H](O)[C@H]8O)CO)[C@H](O)[C@H]7NC(=O)C)CO)CO",  
        # The long SMILES above corresponds to a glycosidic oligomer. In many cases it may miss our criteria
        # because one sugar ring might be fused/linked heavily. (a false negative)
    ]
    for s in test_smiles:
        res, reason = is_amino_sugar(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")