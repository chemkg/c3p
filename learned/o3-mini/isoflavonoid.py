"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: isoflavonoid
Definition: Any 1‐benzopyran with an aryl substituent at position 3.
The function attempts to find a fused bicyclic system where a 6‐membered aromatic benzene ring
is fused (sharing 2 atoms) with a 6‐membered heterocycle that contains an oxygen atom --
and then checks that at least one of the “free” (non‐shared) atoms in the heterocycle bears an aryl substituent.
Note: The process is heuristic and may not capture all isoflavonoids accurately.
"""

from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined (for this purpose) as a molecule containing a 1-benzopyran core --
    that is, a fused bicyclic system where a 6-membered aromatic benzene ring is fused with a 6-membered heterocycle
    containing one oxygen atom -- with an aryl substituent (an extra aromatic ring) attached at the position corresponding to C3.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an isoflavonoid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in molecule"
    
    # Helper functions for ring properties
    def is_benzene(ring):
        # Checks if all atoms are carbon, aromatic and the ring has 6 atoms.
        if len(ring) != 6:
            return False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                return False
        return True

    def has_oxygen(ring):
        # Checks if at least one atom in the ring is oxygen
        return any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)

    # Look for two fused rings (sharing exactly 2 atoms)
    candidate_found = False
    for i in range(len(ring_info)):
        for j in range(i+1, len(ring_info)):
            ring1 = ring_info[i]
            ring2 = ring_info[j]
            shared = set(ring1).intersection(ring2)
            if len(shared) != 2:
                continue  # Require fusion via 2 atoms
            # Require both rings have 6 atoms (typical of benzene and pyran)
            if len(ring1) != 6 or len(ring2) != 6:
                continue

            # Identify a benzopyran candidate: one ring should be a benzene; the other should contain an oxygen.
            if is_benzene(ring1) and has_oxygen(ring2):
                benzene_ring = set(ring1)
                heterocycle = set(ring2)
            elif is_benzene(ring2) and has_oxygen(ring1):
                benzene_ring = set(ring2)
                heterocycle = set(ring1)
            else:
                continue

            # At this point we have a benzopyran candidate.
            # Determine the “free” atoms of the heterocycle (not part of the fusion)
            fusion_atoms = shared
            free_atoms = [atom_idx for atom_idx in heterocycle if atom_idx not in fusion_atoms]

            # Now look for an aryl substituent attached to one of the free atoms.
            # The idea is that one of these free atoms (likely position 3 of the pyran) should have a neighbor
            # that is part of an aromatic ring but not in the heterocycle.
            for atom_idx in free_atoms:
                atom = mol.GetAtomWithIdx(atom_idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in heterocycle:
                        continue  # Skip atoms in the fused core
                    if nbr.GetIsAromatic():
                        # Optionally, we could further check that this neighbor is part of a ring (e.g., a phenyl)
                        # For now, we use aromaticity as a surrogate for "aryl substituent"
                        return True, ("Found benzopyran core (fused ring system with benzene and a pyran-like ring) "
                                      "with an aryl substituent attached at a free position (likely C3).")
    return False, "No benzopyran core with an aryl substituent at position 3 was detected"

# Example usage (you may remove or comment out these lines when integrating):
if __name__ == "__main__":
    # Example: daidzein (an isoflavonoid)
    smiles_examples = [
        "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",  # daidzein
        "CC(=O)O",  # acetyl group (negative control)
    ]
    for smi in smiles_examples:
        is_iso, reason = is_isoflavonoid(smi)
        print(f"SMILES: {smi}\nClassified as isoflavonoid: {is_iso}\nReason: {reason}\n")