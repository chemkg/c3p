"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: CHEBI:3'-hydroxyflavanones
Definition: Any hydroxyflavanone with a hydroxy substituent at position 3'
of the phenyl (B) ring attached to the flavanone core (2-phenylchroman-4-one).
Heuristic:
  1. Find candidate benzene rings (aromatic, exactly 6 atoms).
  2. For atoms in such a benzene ring, check for neighbors not in the ring.
     These neighbors are potential attachment points to the flavanone core.
  3. For each neighbor, look for a ring (other than the candidate benzene ring)
     that is six-membered and non-aromatic. Within that candidate flavanone core ring,
     require exactly one ring oxygen and at least one carbon atom that is bonded to a double-bonded oxygen (ketone).
  4. For the candidate benzene (B) ring, find the atoms meta (two bonds away) 
     from the attachment atom and check that at least one of these meta atoms carries an -OH group.
Notes:
  - These additional criteria help reduce false positives.
"""
from rdkit import Chem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone.
    
    A 3'-hydroxyflavanone is defined as a flavanone (2-phenylchroman-4-one)
    that possesses a hydroxy (-OH) substituent at the meta (3') position relative 
    to the attachment on the B ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a 3'-hydroxyflavanone, False otherwise.
        str: Reason for the classification.
    """
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Helper: Returns True if the given atom has a hydroxyl (-OH) substituent.
    def has_OH(atom):
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # Oxygen
                # Either explicit hydrogen(s) or at least one implicit hydrogen
                if nbr.GetTotalNumHs() > 0:
                    return True
        return False

    # Helper: Given a ring (tuple of atom indices) and attachment index within the ring,
    # return the meta positions (2 bonds away in the ring).
    def meta_indices(ring, att_idx):
        meta_set = set()
        ring_size = len(ring)
        # Locate the position(s) where att_idx occurs (should be one normally)
        positions = [pos for pos, idx in enumerate(ring) if idx == att_idx]
        for pos in positions:
            meta_set.add(ring[(pos + 2) % ring_size])
            meta_set.add(ring[(pos - 2) % ring_size])
        return meta_set

    # Loop over all rings to find candidate aromatic B rings: exactly 6 atoms and fully aromatic.
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue  # not a benzene ring

        # For each atom (attachment candidate) in the benzene ring:
        for att_idx in ring:
            att_atom = mol.GetAtomWithIdx(att_idx)
            # Look at neighbors outside the ring
            neighbor_atoms = [nbr for nbr in att_atom.GetNeighbors() if nbr.GetIdx() not in ring]
            if not neighbor_atoms:
                continue

            # For each neighboring atom, search for a candidate flavanone core ring.
            for nbr in neighbor_atoms:
                candidate_core_found = False
                # Look through rings that include this neighbor but are not the aromatic ring we started with.
                for core_ring in atom_rings:
                    if att_idx in core_ring:
                        # Skip if this ring is the candidate B ring itself.
                        if set(core_ring) == set(ring):
                            continue
                        # For flavanone core, we require a 6-membered ring that is NOT fully aromatic.
                        if len(core_ring) != 6:
                            continue
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in core_ring):
                            continue
                        # Ensure the candidate core ring contains the neighbor atom.
                        if nbr.GetIdx() not in core_ring:
                            continue

                        # Check that in the candidate core ring there is exactly one ring oxygen.
                        oxygen_count = sum(1 for i in core_ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 8)
                        if oxygen_count != 1:
                            continue

                        # Check that at least one carbon in the candidate core ring is double-bonded to an oxygen (ketone).
                        carbonyl_found = False
                        for i in core_ring:
                            atom_core = mol.GetAtomWithIdx(i)
                            if atom_core.GetAtomicNum() == 6:
                                for nbr2 in atom_core.GetNeighbors():
                                    # Get the bond between the carbon and its neighbor.
                                    bond = mol.GetBondBetweenAtoms(atom_core.GetIdx(), nbr2.GetIdx())
                                    if (nbr2.GetAtomicNum() == 8 and 
                                        bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE):
                                        carbonyl_found = True
                                        break
                            if carbonyl_found:
                                break

                        if not carbonyl_found:
                            continue

                        # If we reach here, then this core ring fits our flavanone core criteria.
                        candidate_core_found = True
                        break  # No need to search further for a candidate core ring.

                if not candidate_core_found:
                    continue

                # Once candidate core is found, use the current aromatic ring as the B ring.
                # Find meta positions (2 bonds away in the B ring) from att_idx.
                meta_idxs = meta_indices(ring, att_idx)
                found_meta_OH = False
                for midx in meta_idxs:
                    meta_atom = mol.GetAtomWithIdx(midx)
                    if has_OH(meta_atom):
                        found_meta_OH = True
                        break

                if found_meta_OH:
                    reason = ("Molecule contains a flavanone core "
                              "(a non‐aromatic six‐membered ring with one ring oxygen and a ketone) "
                              "fused to an aromatic B ring that has an -OH substituent at the meta (3') position.")
                    return True, reason

    return False, "Molecule does not appear to contain a 2-phenylchroman-4-one core with an -OH at the 3' position of the B ring"

# Example usage (for testing):
if __name__ == "__main__":
    # Try a few examples.
    examples = [
        ("Oc1cccc(c1)[C@@H]1CC(=O)c2ccccc2O1", "(2S)-3'-hydroxyflavanone (expected True)"),
        ("O1C(C2=CC(=C(O)C(=C2)CC=C(C)C)CC=C(C)C)CC(=O)C3=C1C=C(O)C=C3", "Abyssinone IV (expected False)")
    ]
    for smi, name in examples:
        result, explanation = is_3__hydroxyflavanones(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nExplanation: {explanation}\n")