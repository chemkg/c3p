"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: Nitrohydrocarbon
Definition: A C-nitro compound that is a hydrocarbon in which one or more of the hydrogens
has been replaced by nitro groups.
Heuristic:
  1. The molecule must contain at least one nitro group ([N+](=O)[O-]) that is attached to a carbon.
  2. For each nitro group the attached carbon (found via neighbor connectivity) should show evidence
     of having at least one hydrogen (after explicit H addition) so that the nitro is a substituent.
  3. Removing (i.e. “masking”) all nitro group atoms, the remainder (the backbone) should be predominantly (>=90%)
     carbon atoms.
  4. If the backbone is aromatic (contains rings), then if more than one ring exists the rings must be fused (share at least one bond)
     rather than being isolated “loose” rings.
  
This heuristic is intended to reject cases where extra heteroatoms (outside the nitro groups) occur
or where the hydrocarbon framework is built out of separate ring systems.
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    
    A nitrohydrocarbon must have at least one nitro group ([N+](=O)[O-]) attached directly to a carbon.
    Furthermore, when we remove all atoms that are part of the nitro groups, the remaining heavy atoms
    (i.e. non-hydrogens) should be (almost) all carbons, and if aromatic rings are present these should
    be part of a single fused system (not, for instance, two aromatic rings connected by a single bond).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the nitro group SMARTS pattern.
    nitro_smarts = Chem.MolFromSmarts("[N+](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_smarts):
        return False, "No nitro group ([N+](=O)[O-]) found in the molecule"
    
    # Check that at least one nitro group is attached to a carbon.
    nitro_matches = mol.GetSubstructMatches(nitro_smarts)
    attached_to_carbon = False
    # For checking if a hydrogen was replaced, we will add explicit hydrogens.
    molH = Chem.AddHs(mol)
    for match in nitro_matches:
        # In each nitro match the atoms are (typically) [N, O, O].
        # Find the nitrogen and then look for its neighbor that is carbon and not in the nitro group.
        nitroN = match[0]  # assume first atom is the nitrogen
        atomN = mol.GetAtomWithIdx(nitroN)
        for nbr in atomN.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in match:
                # Check that in the hydrogen-added molecule, this carbon has at least one H.
                nbrH = molH.GetAtomWithIdx(nbr.GetIdx())
                if nbrH.GetTotalNumHs() > 0:
                    attached_to_carbon = True
                    break
        if attached_to_carbon:
            break
    if not attached_to_carbon:
        return False, "Nitro group(s) not directly attached to a carbon atom that had a hydrogen"

    # Collect indices of all atoms that are part of any nitro group.
    nitro_atom_indices = set()
    for match in nitro_matches:
        nitro_atom_indices.update(match)
    
    # Create a copy of the molecule (backbone) by removing nitro group atoms.
    # We work on an editable copy.
    backbone = Chem.RWMol(mol)
    # Remove in descending order of indices to not mess up coordinates.
    for idx in sorted(nitro_atom_indices, reverse=True):
        try:
            backbone.RemoveAtom(idx)
        except Exception:
            pass
    backbone = backbone.GetMol()
    
    # Get heavy atoms (non-hydrogen) in the backbone.
    heavy_atoms = [atom for atom in backbone.GetAtoms() if atom.GetAtomicNum() != 1]
    if not heavy_atoms:
        return False, "Backbone is empty after removing nitro groups"
    
    # Count how many heavy atoms are carbon.
    n_carbon = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
    hydro_ratio = n_carbon / len(heavy_atoms)
    if hydro_ratio < 0.9:
        return False, f"Backbone has only {hydro_ratio*100:.1f}% carbons (requires at least 90% C); extra heteroatoms present"
    
    # For aromatic backbones: if rings are present, check that rings are fused.
    ri = backbone.GetRingInfo()
    rings = ri.AtomRings()
    if len(rings) > 1:
        # For each ring, determine bonds (as frozensets of atom index pairs)
        ring_bonds = []
        for ring in rings:
            bonds = set()
            ring = list(ring)
            n = len(ring)
            for i in range(n):
                a1 = ring[i]
                a2 = ring[(i+1)%n]
                bond = frozenset((a1, a2))
                bonds.add(bond)
            ring_bonds.append(bonds)
        fused = False
        # if at least one pair of rings shares a bond, we consider them fused
        for i in range(len(ring_bonds)):
            for j in range(i+1, len(ring_bonds)):
                if ring_bonds[i].intersection(ring_bonds[j]):
                    fused = True
                    break
            if fused:
                break
        if not fused:
            return False, "Aromatic backbone rings are not fused (found separate isolated rings)"
    
    # Passed all the checks
    return True, "Molecule is a nitrohydrocarbon: nitro group(s) attached to a (largely) hydrocarbon backbone"

# Uncomment below to do quick tests.
# test_smiles_list = [
#     "[O-][N+](=O)c1ccc2c3ccccc3c2c1",  # 1-nitronaphthalene (true positive expected)
#     "[O-][N+](=O)c1ccc(cc1)-c1ccccc1",   # 4-Nitrobiphenyl (false positive expected)
# ]
# for s in test_smiles_list:
#     result, reason = is_nitrohydrocarbon(s)
#     print(s, result, reason)