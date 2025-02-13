"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: 17beta-hydroxy steroid
Definition: A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.
This routine heuristically checks that the molecule has a fused steroid nucleus (three 6-membered rings and one 5-membered ring)
and that at least one chiral carbon in a 5-membered ring carries an –OH group.
"""

from rdkit import Chem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    This is a heuristic approach that:
      1. Verifies the presence of a steroid-like fused ring system (one five-membered and at least three six-membered rings).
      2. Searches for a chiral carbon (assumed to be at the 17-position) in a five-membered ring that bears an –OH group.
         (We assume that in a conventional steroid depiction, the 17beta-OH will be encoded as a chiral center with a specific tag.)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 17beta-hydroxy steroid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Retrieve ring information: get all rings (each as a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Count rings by size. Steroid nucleus is typically three 6-membered rings and one 5-membered ring.
    count_5 = sum(1 for ring in ring_info if len(ring) == 5)
    count_6 = sum(1 for ring in ring_info if len(ring) == 6)
    if count_5 < 1 or count_6 < 3:
        return False, ("Molecule does not appear to have a typical fused steroid core "
                       "(expected at least one 5-membered and three 6-membered rings)")

    # Search for candidate chiral carbon in a 5-membered ring that carries an –OH group.
    candidate_found = False
    for atom in mol.GetAtoms():
        # Focus on carbon atoms that have chirality information
        if atom.GetAtomicNum() != 6:
            continue
        # Check if the atom is defined as chiral
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue

        # Look for an attached hydroxyl group (an oxygen connected by a SINGLE bond)
        oh_found = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    oh_found = True
                    break
        if not oh_found:
            continue

        # Check if atom is part of at least one 5-membered ring.
        # (The steroid D ring is five-membered where C-17 is usually located.)
        atom_rings_5 = [set(ring) for ring in ring_info if atom.GetIdx() in ring and len(ring) == 5]
        if not atom_rings_5:
            continue

        # (Optional) Increase our confidence by requiring that the 5-membered ring is fused with at least one 6-membered ring.
        fused_with_6 = False
        for ring5 in atom_rings_5:
            for ring in ring_info:
                if len(ring) == 6:
                    if len(ring5.intersection(ring)) >= 2:
                        fused_with_6 = True
                        break
            if fused_with_6:
                break
        if not fused_with_6:
            continue

        # At this point, we found a chiral carbon in a 5-membered ring carrying an –OH and fused with a six-membered ring.
        # We assume this corresponds to the 17beta-hydroxy group.
        candidate_found = True
        break

    if candidate_found:
        return True, "Molecule contains a steroid nucleus with a candidate 17beta-hydroxy group."
    else:
        return False, "No candidate 17beta-hydroxy group was detected in the steroid framework."

# Example usage:
if __name__ == '__main__':
    # Test with one of the provided examples: 17beta-estradiol 
    test_smiles = "[H][C@]12CC[C@]3(C)[C@@H](O)CC[C@@]3([H])[C@]1([H])CCc1cc(O)ccc21"
    result, reason = is_17beta_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)