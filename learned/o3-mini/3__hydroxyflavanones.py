"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: CHEBI:3'-hydroxyflavanones
Definition: Any hydroxyflavanone with a hydroxy substituent at position 3'
of the phenyl (B) ring attached to the flavanone core (2-phenylchroman-4-one)
Heuristic:
  1. Find an aromatic benzene ring (6 atoms) as the candidate B ring.
  2. Identify an atom in that ring that attaches to another (non-aromatic) ring,
     candidate for the flavanone core.
  3. In that attached ring (flavanone core), check for a carbonyl group (C=O) and a ring oxygen.
  4. For the benzene ring, determine the atoms that are meta (two bonds away) from the attachment.
  5. Check that at least one meta atom carries a hydroxyl (-OH) substituent.
Notes:
  - We fixed the error where GetIsInRing was used and replaced it with IsInRing().
  - The code now leverages RDKit functions correctly to search for the target motifs.
"""
from rdkit import Chem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone.

    A 3'-hydroxyflavanone is defined as a flavanone (2-phenylchroman-4-one)
    that possesses a hydroxy (-OH) substituent at the 3' position (meta relative 
    to the attachment on the B ring) of the phenyl (B) ring.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a 3'-hydroxyflavanone, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Helper: checks if an atom has a hydroxyl (-OH) substituent
    def has_OH(atom):
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Check if oxygen has at least one hydrogen (implicit or explicit)
                if nbr.GetTotalNumHs() > 0:
                    return True
        return False

    # Helper: given a ring (as a tuple/list of atom indices) and an attachment index within the ring,
    # compute the meta positions (two bonds away in the ring).
    def meta_indices(ring, att_idx):
        positions = [i for i, idx in enumerate(ring) if idx == att_idx]
        meta_set = set()
        ring_size = len(ring)
        for pos in positions:
            meta_set.add(ring[(pos + 2) % ring_size])
            meta_set.add(ring[(pos - 2) % ring_size])
        return meta_set

    # Look for candidate benzene rings: aromatic, exactly 6 members.
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue

        # For each atom in the benzene ring, see if it attaches to another ring (candidate core)
        for att_idx in ring:
            att_atom = mol.GetAtomWithIdx(att_idx)
            # Consider neighbors outside the benzene ring
            neighbors_outside = [nbr for nbr in att_atom.GetNeighbors() if nbr.GetIdx() not in ring]
            if not neighbors_outside:
                continue

            for nbr in neighbors_outside:
                # Look at all rings that contain the neighbor (candidate flavanone core ring)
                candidate_core = False
                for core_ring in atom_rings:
                    # The neighbor must be in the candidate core ring and the ring must not be fully aromatic.
                    if nbr.GetIdx() not in core_ring:
                        continue
                    if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in core_ring):
                        continue

                    # Check for the carbonyl (C=O) group and a ring oxygen in that ring.
                    carbonyl_found = False
                    ring_oxygen_found = False
                    for idx in core_ring:
                        atom_core = mol.GetAtomWithIdx(idx)
                        # Check if the atom is oxygen and part of the ring.
                        if atom_core.GetAtomicNum() == 8 and atom_core.IsInRing():
                            ring_oxygen_found = True
                        # For a carbon, look for a double-bonded oxygen.
                        if atom_core.GetAtomicNum() == 6:
                            for nbr2 in atom_core.GetNeighbors():
                                bond = mol.GetBondBetweenAtoms(atom_core.GetIdx(), nbr2.GetIdx())
                                if nbr2.GetAtomicNum() == 8 and bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                                    carbonyl_found = True
                    if carbonyl_found and ring_oxygen_found:
                        candidate_core = True
                        break

                if not candidate_core:
                    continue

                # Now treat the candidate benzene ring as the B ring.
                # Get the meta positions relative to the attachment atom.
                meta_idxs = meta_indices(ring, att_idx)
                found_meta_OH = False
                for midx in meta_idxs:
                    meta_atom = mol.GetAtomWithIdx(midx)
                    if has_OH(meta_atom):
                        found_meta_OH = True
                        break

                if found_meta_OH:
                    reason = ("Molecule contains a flavanone core (chromanone with a carbonyl and ring oxygen) "
                              "and a B ring (aromatic 6-membered ring) with a hydroxyl substituent at the meta (3') position.")
                    return True, reason

    return False, "Molecule does not appear to contain a 2-phenylchroman-4-one core with an -OH at the 3' position of the B ring"

# Example usage: Feel free to test with provided SMILES examples.
if __name__ == "__main__":
    # (2S)-3'-hydroxyflavanone (example)
    test_smiles = "Oc1cccc(c1)[C@@H]1CC(=O)c2ccccc2O1"
    result, explanation = is_3__hydroxyflavanones(test_smiles)
    print(result, explanation)