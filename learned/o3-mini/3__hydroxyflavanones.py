"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: CHEBI:3'-hydroxyflavanones
Definition: Any hydroxyflavanone with a hydroxy substituent at position 3'
of the phenyl (B) ring attached to the flavanone core.
Note: The detection is heuristic. It first checks for a flavanone (2-phenylchroman-4-one)
core by a SMARTS pattern and then examines the phenyl ring attached at C2 to see if it
contains a hydroxyl group in the meta position (2 bonds away from the attachment).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone.
    
    A 3'-hydroxyflavanone is defined as a flavanone (a 2-phenylchroman-4-one)
    that possesses a hydroxy substituent at the 3' position of the phenyl (B) ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, check for a flavanone (chroman-4-one) core.
    # The following SMARTS captures the basic 2-phenylchroman-4-one skeleton.
    # It looks for a six‐membered ring containing a carbonyl and a ring‐oxygen.
    flavanone_core_smarts = "C1CC(=O)Oc2ccccc12"
    core_pattern = Chem.MolFromSmarts(flavanone_core_smarts)
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Molecule does not contain a flavanone (2-phenylchroman-4-one) core"

    # Next, we try to locate the B ring (the non-fused phenyl ring attached to C2).
    # We do this by iterating over bonds from sp3 carbons in the fused ring system
    # (the chromanone core) to an aromatic carbon that is not in that same fused system.
    matches = mol.GetSubstructMatches(core_pattern)
    # We expect at least one match for the flavanone core; take the first.
    core_atoms = set(matches[0])
    # Find candidate bonds: an sp3 carbon (likely the C2 chiral center) in the core that
    # is connected to an aromatic atom that is NOT part of the core.
    candidate_attachment = None
    for atom in mol.GetAtoms():
        # look for an sp3 carbon that is in the core (part of the chromanone ring)
        # (we use hybridization to require it is sp3)
        if atom.GetAtomicNum() == 6 and atom.GetHybridization().name == "SP3" and atom.GetIdx() in core_atoms:
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in core_atoms and nbr.GetIsAromatic():
                    candidate_attachment = (atom, nbr)
                    break
            if candidate_attachment:
                break

    if candidate_attachment is None:
        return False, "Flavanone core found, but could not locate the attached B ring"

    # candidate_attachment[1] is the aromatic atom (in the B ring) attached to the core.
    attachment_atom = candidate_attachment[1]

    # Now, find a 6-membered aromatic ring (benzene) that contains this attachment atom.
    ring_found = False
    ring_info = mol.GetRingInfo()
    b_ring = None
    for ring in ring_info.AtomRings():
        # Look for rings of size 6 in which every atom is aromatic.
        if len(ring) == 6:
            if attachment_atom.GetIdx() in ring:
                # Verify all atoms in the ring are aromatic.
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    b_ring = ring
                    ring_found = True
                    break

    if not ring_found or b_ring is None:
        return False, "B ring (aromatic 6-membered ring connected to core) not found"

    # In the benzene ring, the attachment_atom is at the 1' position.
    # For a hydroxyl to be at the 3' position, it should be attached
    # to an atom that is two bonds away from the attachment atom along the ring.
    # We now examine the atoms in b_ring for an -OH substituent.
    found_meta_OH = False
    for idx in b_ring:
        # Skip the attachment point.
        if idx == attachment_atom.GetIdx():
            continue
        candidate = mol.GetAtomWithIdx(idx)
        # Check if candidate has a hydroxyl group directly attached.
        # Look at neighbors: an oxygen with at least one hydrogen (implicit or explicit).
        for nbr in candidate.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Check if it looks like a hydroxyl (not a carbonyl or ether).
                # We require that the oxygen has at least one hydrogen.
                if nbr.GetTotalNumHs() > 0:
                    # Now check if candidate is meta to the attachment_atom on the ring.
                    # Compute the shortest path between attachment_atom and candidate.
                    path = rdmolops.GetShortestPath(mol, attachment_atom.GetIdx(), candidate.GetIdx())
                    # In a benzene ring the meta relationship gives a path length of 3 atoms (2 bonds).
                    if len(path) == 3:
                        found_meta_OH = True
                        break
        if found_meta_OH:
            break

    if not found_meta_OH:
        return False, "Flavanone core found but no hydroxy substituent at the 3' (meta) position of the B ring"

    return True, "Molecule contains a flavanone core with a hydroxy group meta (3'-position) on the B ring"

# Example usage (you can remove or comment out these lines when integrating):
if __name__ == "__main__":
    test_smiles = "Oc1cccc(c1)[C@@H]1CC(=O)c2ccccc2O1"  # (2S)-3'-hydroxyflavanone
    result, reason = is_3__hydroxyflavanones(test_smiles)
    print(result, reason)