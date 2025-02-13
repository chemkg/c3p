"""
Classifies: CHEBI:28963 amino sugar
"""
#!/usr/bin/env python
"""
Classifies amino sugars.
Definition: An amino sugar is a sugar molecule in which one or more alcoholic hydroxyl groups 
            attached to the sugar‐ring have been replaced by an unsubstituted amino group.

Algorithm:
  1. Parse the SMILES string.
  2. Identify candidate rings by:
       • Considering only rings of size 5 or 6.
       • Requiring that the ring is non‐aromatic, contains exactly one oxygen and the other atoms are carbons.
       • Rejecting rings in which any atom is present in more than one ring (fused rings).
  3. For each candidate ring, look at substituents attached to ring carbons (ignoring the ring oxygen):
       • Only consider small one‐atom substituents (either single atom neighbors or those with only H’s apart from the ring bond).
       • Count a substituent as a hydroxyl if the neighbor is oxygen with at least one explicit hydrogen.
       • Count as “free amino” only if the neighbor is a nitrogen with exactly two explicit hydrogens and no double‐bond to oxygen.
       • A candidate ring must have at least three polar substituents (–OH or free –NH₂) and at least one of them must be a free –NH₂.
  4. Return True if at least one candidate ring meets the criteria.
  
If no candidate ring is found with the required free amino substituent, return False.
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
    
    # Update property cache to ensure rings are computed properly.
    mol.UpdatePropertyCache(strict=False)
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # tuple of tuples of atom indices

    candidate_rings = []  # each candidate as a set of indices
    # Loop over each ring from RDKit
    for ring in all_rings:
        # Accept only 5- or 6-membered rings
        if len(ring) not in [5, 6]:
            continue

        # Get the atoms for the ring
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # We require the ring to be non-aromatic.
        if any(atom.GetIsAromatic() for atom in ring_atoms):
            continue

        # Count oxygen atoms and ensure all nonoxygen atoms are carbons.
        oxygen_count = 0
        valid_atoms = True
        for atom in ring_atoms:
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            elif atom.GetAtomicNum() != 6:
                valid_atoms = False
                break
        if not valid_atoms or oxygen_count != 1:
            continue

        # Check that none of the atoms in this ring are shared with another ring (i.e. a fused ring)
        fused = False
        for idx in ring:
            # Count how many rings contain this atom.
            count = sum(1 for r in all_rings if idx in r)
            if count > 1:
                fused = True
                break
        if fused:
            continue
        
        candidate_rings.append(set(ring))
        
    if not candidate_rings:
        return False, "No isolated sugar-like ring (5- or 6-membered, non-aromatic with exactly one oxygen) detected."
    
    # Helper: Check if an atom qualifies as a hydroxyl (-OH) substituent.
    def is_hydroxyl(neighbor):
        if neighbor.GetAtomicNum() != 8:
            return False
        if neighbor.IsInRing():
            return False
        # Check explicit hydrogens on the oxygen (at least one required)
        if neighbor.GetTotalNumHs(includeNeighbors=False) < 1:
            return False
        return True

    # Helper: Check if an atom qualifies as a free amino (-NH2) substituent.
    def is_free_amino(neighbor):
        if neighbor.GetAtomicNum() != 7:
            return False
        if neighbor.IsInRing():
            return False
        # Require exactly 2 explicit hydrogens
        if neighbor.GetTotalNumHs(includeNeighbors=False) != 2:
            return False
        # Check that the nitrogen is not directly involved in a C=O double bond (which would suggest an amide)
        for bond in neighbor.GetBonds():
            # RDKit bond type for a double bond is accessed via GetBondType()
            # We check if the bond is double and the other atom is oxygen.
            if bond.GetBondTypeAsDouble() == 2.0:
                other = bond.GetOtherAtom(neighbor)
                if other.GetAtomicNum() == 8:
                    return False
        return True

    # Now check each candidate ring for the correct substituents.
    for ring_set in candidate_rings:
        polar_substituent_count = 0
        found_free_amino = False

        # Consider only ring atoms that are carbons (skip the ring oxygen)
        for idx in ring_set:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Check each neighbor of the ring carbon that is not part of the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                # Only consider small substituents, meaning that the neighbor should either have degree 1
                # or, if its degree is larger, all bonds (except the bond attached to our ring atom) are to hydrogens.
                if nbr.GetDegree() > 1:
                    heavy_neighbors = 0
                    for nn in nbr.GetNeighbors():
                        if nn.GetIdx() == atom.GetIdx():
                            continue
                        if nn.GetAtomicNum() != 1:
                            heavy_neighbors += 1
                    if heavy_neighbors > 0:
                        continue

                # Test if the substituent qualifies as hydroxyl or free amino.
                if is_hydroxyl(nbr):
                    polar_substituent_count += 1
                elif is_free_amino(nbr):
                    polar_substituent_count += 1
                    found_free_amino = True
        # A candidate ring must have at least three polar substituents and at least one free amino.
        if found_free_amino and polar_substituent_count >= 3:
            return True, ("Molecule contains an isolated sugar-like ring (5- or 6-membered with one oxygen) " +
                          "decorated with at least three small polar substituents (including a free -NH2) " +
                          "consistent with an amino sugar structure.")
    
    return False, "No candidate sugar ring with a free -NH2 replacing an -OH group was detected."


# Example usage:
if __name__ == "__main__":
    test_smiles = [
        # 3-amino-3-deoxy-D-glucopyranose; expected True.
        "N[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O",
        # A complex sugar-like oligomer (likely not matching our strict criteria)
        "O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H](O)[C@H]1O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)CO[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O[C@@H]7O[C@@H]([C@@H](O[C@@H]8O[C@@H]([C@H](O)[C@H](O)[C@H]8O)CO)[C@H](O)[C@H]7NC(=O)C)CO)CO"
    ]
    for s in test_smiles:
        res, reason = is_amino_sugar(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")