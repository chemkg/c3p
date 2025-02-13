"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: CHEBI:3'-hydroxyflavanones
Definition: Any hydroxyflavanone with a hydroxy substituent at position 3'
of the phenyl (B) ring attached to the flavanone core.
Heuristic:
  1. Identify an aromatic benzene ring (a candidate B ring) with 6 atoms.
  2. Find an atom in that ring that attaches to another (non‐aromatic) ring.
  3. In that attached ring, verify that there is a carbonyl group (C=O) and a ring oxygen,
     which are key features of the chromanone (flavanone) core.
  4. In the benzene ring, determine the atoms that are meta (two bonds away) from the attachment.
  5. Check that at least one of those meta atoms carries a hydroxyl group.
  
This two‐step approach is more tolerant of substituents while aiming to capture the 2-phenylchroman-4-one core.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone.
    
    A 3'-hydroxyflavanone is defined as a flavanone (2-phenylchroman-4-one)
    that possesses a hydroxy (-OH) substituent at the 3' position (meta relative 
    to the attachment on the B ring) of the phenyl (B) ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    candidate_found = False

    # Helper function: check if an atom has a hydroxyl (-OH) substituent.
    def has_OH(atom):
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Check if this oxygen has at least one hydrogen
                # (RDKit may store this as implicit or explicit hydrogens)
                if nbr.GetTotalNumHs() > 0:
                    return True
        return False

    # Helper function: determine meta indices in a benzene ring.
    # We assume the ring is a cyclic list. We try every cyclic permutation.
    def meta_indices(ring, att_idx):
        # Find the position(s) within ring that equal att_idx.
        positions = [i for i, idx in enumerate(ring) if idx == att_idx]
        meta_set = set()
        ring_size = len(ring)
        for pos in positions:
            # In a cyclic order, atoms two positions away are meta.
            meta_set.add(ring[(pos + 2) % ring_size])
            meta_set.add(ring[(pos - 2) % ring_size])
        return meta_set

    # We now iterate over all rings that are aromatic and have exactly 6 atoms (candidate benzene rings)
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        # Verify that every atom in the ring is aromatic.
        if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue

        # For each atom in the benzene ring, find a neighbor that is not in the ring.
        # This neighbor is a candidate for the attachment to the flavanone core.
        for att_idx in ring:
            att_atom = mol.GetAtomWithIdx(att_idx)
            neighbors_outside = [nbr for nbr in att_atom.GetNeighbors() if nbr.GetIdx() not in ring]
            if not neighbors_outside:
                continue

            # For each neighbor outside the ring, check if it belongs to a candidate core ring.
            for nbr in neighbors_outside:
                # We want the neighbor to be in a ring that is not fully aromatic.
                # (Often the chromanone core is a saturated ring or partially unsaturated.)
                nbr_rings = nbr.GetOwningMol().GetRingInfo().AtomRings()
                # We can check each ring of which nbr is a part:
                core_candidate = False
                for r in nbr_rings:
                    if nbr.GetIdx() not in r:
                        continue
                    # If at least one atom in the ring is non‐aromatic, consider it.
                    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r):
                        # Check if a carbonyl group is present in this ring.
                        # Look for a carbon (atomic num 6) with a double-bonded oxygen 
                        # (i.e. a neighbor oxygen with bond type DOUBLE).
                        carbonyl_found = False
                        ring_has_oxygen = False
                        for i in r:
                            atom_core = mol.GetAtomWithIdx(i)
                            if atom_core.GetAtomicNum() == 8 and atom_core.GetIsInRing():
                                ring_has_oxygen = True
                            if atom_core.GetAtomicNum() == 6:
                                for nbr2 in atom_core.GetNeighbors():
                                    # Check for C=O bond by looking at bond order.
                                    bond = mol.GetBondBetweenAtoms(atom_core.GetIdx(), nbr2.GetIdx())
                                    if nbr2.GetAtomicNum() == 8 and bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                                        carbonyl_found = True
                        if carbonyl_found and ring_has_oxygen:
                            core_candidate = True
                            break

                if not core_candidate:
                    continue

                # We now treat the current benzene ring as the B ring and the att_atom as the attachment.
                # Determine which atoms in the ring are meta relative to the attachment.
                meta_idxs = meta_indices(ring, att_idx)
                # Check if any meta atom in the ring carries an -OH group.
                found_meta_OH = False
                for midx in meta_idxs:
                    meta_atom = mol.GetAtomWithIdx(midx)
                    if has_OH(meta_atom):
                        found_meta_OH = True
                        break

                if found_meta_OH:
                    reason = ("Molecule contains a candidate flavanone core (a non‐aromatic ring with a carbonyl and ring‐oxygen) "
                              "and a B ring (aromatic 6‐membered ring) attached at an atom whose meta position carries a hydroxyl group")
                    return True, reason
    # If we finish looping without finding a matching candidate:
    return False, "Molecule does not appear to contain a 2-phenylchroman-4-one core with a hydroxy substituent at the 3' (meta) position of the B ring"

# Example usage (feel free to test with provided SMILES):
if __name__ == "__main__":
    # (2S)-3'-hydroxyflavanone (example)
    test_smiles = "Oc1cccc(c1)[C@@H]1CC(=O)c2ccccc2O1"
    result, reason = is_3__hydroxyflavanones(test_smiles)
    print(result, reason)