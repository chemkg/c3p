"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: Catechin class – members of the hydroxyflavan group that contain a flavan-3-ol (catechin) scaffold.
Heuristic (improved):
  • Convert the SMILES to a molecule and add explicit hydrogens so that –OH groups are detectable.
  • Iterate over all 6‐membered rings and look for a candidate “C ring” that:
      - Contains exactly one ring oxygen.
      - Has at least three atoms with sp3 hybridization (to avoid fully aromatic rings).
  • Check that at least one sp3 carbon (not the ring oxygen) has a hydroxyl (-OH) substituent.
  • Verify that the candidate ring is fused with at least one aromatic ring (“A ring”) 
    by sharing at least 2 atoms.
  • Check that at least one candidate ring atom is bonded to an aromatic ring (“B ring”)
    not part of the fused aromatic system.
Because of the diversity in substitution patterns, this algorithm is heuristic.
"""

from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule belongs to the catechin (flavan-3-ol) class based on its SMILES.
    The heuristic searches for a candidate 6-membered (C) ring with one oxygen and several sp3 atoms,
    a hydroxyl substituent on a candidate carbon, fused to an aromatic ring (A ring) and having an
    additional aromatic substituent (B ring).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a catechin scaffold, False otherwise.
        str: Reason for the classification decision.
    """
    # Try to parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    # Add explicit hydrogens to detect hydroxyl groups.
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    
    # Helper: determine if a ring (list/tuple of atom indices) is entirely aromatic.
    def is_aromatic_ring(ring):
        return all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
    
    # Iterate over candidate rings.
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # we only want 6-membered rings
        candidate_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        
        # Criterion 1: Candidate ring should contain exactly one oxygen.
        oxygen_count = sum(1 for atom in candidate_atoms if atom.GetAtomicNum() == 8)
        if oxygen_count != 1:
            continue
        
        # Criterion 2: At least three atoms in the candidate ring should be sp3.
        sp3_count = sum(1 for atom in candidate_atoms if atom.GetHybridization() == Chem.HybridizationType.SP3)
        if sp3_count < 3:
            continue
        
        # Criterion 3: At least one carbon (atomic number 6) in the candidate ring should have an -OH substituent.
        has_oh = False
        for atom in candidate_atoms:
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                # Only consider neighbors that are not in the candidate ring.
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() >= 1:
                    has_oh = True
                    break
            if has_oh:
                break
        if not has_oh:
            # If no hydroxyl is found on any candidate carbon, skip to the next candidate.
            continue
        
        # Criterion 4: Find a fused aromatic ring (A ring).
        fused_aromatic = None
        for other_ring in atom_rings:
            if set(other_ring) == set(ring):
                continue
            common_atoms = set(ring).intersection(other_ring)
            if len(common_atoms) >= 2 and is_aromatic_ring(other_ring):
                fused_aromatic = set(other_ring)
                break
        if fused_aromatic is None:
            continue

        # Criterion 5: Look for at least one substituent aromatic ring (B ring) attached to a candidate ring atom
        # that is not already part of the fused aromatic system.
        b_ring_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring or (fused_aromatic and nbr_idx in fused_aromatic):
                    continue
                # Check every ring that neighbor participates in.
                for other_ring in atom_rings:
                    if nbr_idx in other_ring:
                        # if this ring is not identical to the fused aromatic ring and is aromatic, count it as the B-ring.
                        if fused_aromatic and set(other_ring) == fused_aromatic:
                            continue
                        if is_aromatic_ring(other_ring):
                            b_ring_found = True
                            break
                if b_ring_found:
                    break
            if b_ring_found:
                break
        if not b_ring_found:
            continue

        # If we reached here, then our candidate ring fulfilled our heuristic criteria
        # for having a flavan-3-ol (catechin) scaffold.
        return True, ("Molecule contains a catechin (flavan-3-ol) scaffold: a 6-membered oxygen-containing "
                      "ring (with predominantly sp3 atoms and a hydroxyl substituent), fused with an aromatic "
                      "ring and bearing an additional aromatic substituent.")

    # If no candidate ring passes the filters, return failure.
    return False, "Molecule does not appear to contain the required flavan-3-ol (catechin) scaffold."

# Example usage:
if __name__ == "__main__":
    # Test with one example: (-)-catechin 
    test_smiles = "O[C@@H]1Cc2c(O)cc(O)cc2O[C@H]1c1ccc(O)c(O)c1"
    decision, reason = is_catechin(test_smiles)
    print("Is catechin?", decision)
    print("Reason:", reason)