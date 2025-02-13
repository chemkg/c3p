"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: Nitrohydrocarbon
Definition: A C-nitro compound that is a hydrocarbon in which one or more hydrogens has been replaced by nitro groups.
Heuristic:
  1. The molecule must contain at least one nitro group ([N+](=O)[O-]) attached directly to a carbon atom.
  2. Remove all atoms belonging to nitro groups; the remaining heavy atoms must be almost exclusively (>=90%) carbon.
     However, a very limited number of oxygen atoms are allowed—for example, if the oxygen is part of a three‐membered ring (epoxide) that may occur in otherwise all‐carbon fused ring systems.
  3. For aromatic multi‐ring systems, if more than one ring exists, the rings should be fused.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    
    A nitrohydrocarbon is defined as a compound containing one or more nitro groups ([N+](=O)[O-])
    that are attached directly to a carbon, and whose (non-nitro) backbone is predominantly a hydrocarbon.
    
    The heuristic in this implementation is as follows:
    
      1. Check for the presence of at least one nitro group.
      2. Verify that at least one nitro group is attached to a carbon atom.
      3. Remove all atoms that are part of any nitro group.
      4. Examine the remaining (backbone) heavy atoms. They are allowed only if:
            - They are carbon atoms, or
            - They are oxygen atoms that look “epoxide-like” (i.e. in a ring and with exactly two neighbors).
      5. The backbone must contain at least 90% carbon atoms (by count of allowed heavy atoms).
      6. If the backbone is aromatic and contains multiple rings, require that the aromatic rings are fused (i.e. share at least one bond).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a nitrohydrocarbon, False otherwise.
        str: A reason describing the outcome.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the nitro group SMARTS pattern.
    nitro_smarts = Chem.MolFromSmarts("[N+](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_smarts):
        return False, "No nitro group ([N+](=O)[O-]) found in the molecule"
    
    # Check that at least one nitro group is attached directly to a carbon atom.
    nitro_matches = mol.GetSubstructMatches(nitro_smarts)
    attached_to_carbon = False
    for match in nitro_matches:
        # Assume the first atom in the nitro match is the nitrogen.
        nitroN_idx = match[0]
        atomN = mol.GetAtomWithIdx(nitroN_idx)
        # Look at neighbors of the nitro nitrogen.
        for neighbor in atomN.GetNeighbors():
            # If a neighbor is carbon and is not part of this found nitro group, mark valid.
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in match:
                attached_to_carbon = True
                break
        if attached_to_carbon:
            break
    if not attached_to_carbon:
        return False, "Nitro group(s) not directly attached to a carbon atom"
    
    # Gather indices of all atoms that are part of any nitro group.
    nitro_atom_indices = set()
    for match in nitro_matches:
        nitro_atom_indices.update(match)
    
    # Create a copy of the molecule and remove nitro group atoms to obtain the backbone.
    backbone = Chem.RWMol(mol)
    # Remove atoms in descending order so indices do not shift.
    for idx in sorted(nitro_atom_indices, reverse=True):
        try:
            backbone.RemoveAtom(idx)
        except Exception:
            pass
    backbone = backbone.GetMol()
    
    # Examine the backbone heavy atoms (exclude hydrogens).
    backbone_atoms = [atom for atom in backbone.GetAtoms() if atom.GetAtomicNum() != 1]
    if not backbone_atoms:
        return False, "Backbone is empty after removing nitro groups"
    
    # Allowed atoms: carbon always; oxygen allowed only if it appears to be in an epoxide-like environment.
    allowed_atom_types = {6}  # carbon atomic number is 6
    non_carbon_atoms = 0
    carbon_count = 0
    for atom in backbone_atoms:
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:
            carbon_count += 1
        elif atomic_num == 8:
            # Check if oxygen is allowed (for example, in an epoxide-like ring):
            # Typically, an oxygen in an epoxide is part of a ring and has exactly two heavy-atom neighbors.
            if atom.IsInRing() and len([nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1])==2:
                non_carbon_atoms += 1
            else:
                return False, "Backbone contains an oxygen atom not in an allowed (e.g. epoxide) environment"
        else:
            # Any other heteroatom (e.g. N, S, etc.) is not allowed
            return False, f"Backbone contains disallowed heteroatom with atomic number {atomic_num}"
    
    total_allowed = len(backbone_atoms)
    hydro_ratio = carbon_count / total_allowed
    if hydro_ratio < 0.9:
        return False, f"Backbone has only {hydro_ratio*100:.1f}% carbons (requires at least 90%); extra heteroatoms present"
    
    # For aromatic systems with multiple rings, verify that the rings are fused
    ri = backbone.GetRingInfo()
    rings = ri.AtomRings()
    if len(rings) > 1:
        # For each ring, get the set of bonds that form the ring.
        ring_bonds = []
        for ring in rings:
            bonds = set()
            ring = list(ring)
            n = len(ring)
            for i in range(n):
                a1 = ring[i]
                a2 = ring[(i + 1) % n]
                bond = frozenset((a1, a2))
                bonds.add(bond)
            ring_bonds.append(bonds)
        # Check that at least one pair of rings share at least one bond.
        fused = False
        for i in range(len(ring_bonds)):
            for j in range(i + 1, len(ring_bonds)):
                if ring_bonds[i].intersection(ring_bonds[j]):
                    fused = True
                    break
            if fused:
                break
        if not fused:
            return False, "Aromatic backbone rings are not fused (found separate isolated rings)"
    
    return True, "Molecule is a nitrohydrocarbon: nitro group(s) are attached directly to a carbon and the backbone is predominantly hydrocarbon"

# Uncomment and run the following block to perform quick tests:
# test_smiles = [
#     "[O-][N+](=O)c1cccc2C3OC3C=Cc12",  # 1-Nitronaphthalene-5,6-oxide (expected: True)
#     "[O-][N+](=O)c1ccc-2c(c1)-c1ccc([N+]([O-])=O)c3cccc-2c13",  # 3,9-Dinitrofluoranthene (expected: True)
#     "CC(C)[N+]([O-])=O",  # 2-nitropropane (expected: True)
#     "[O-][N+](=O)c1ccc(cc1)-c1ccccc1",  # 4-Nitrobiphenyl (expected: False, extra aromatic system not exclusively hydrocarbon)
# ]
#
# for s in test_smiles:
#     valid, reason = is_nitrohydrocarbon(s)
#     print(f"SMILES: {s}\nResult: {valid}\nReason: {reason}\n")