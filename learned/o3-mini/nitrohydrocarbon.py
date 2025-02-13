"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: Nitrohydrocarbon
Definition: A C-nitro compound that is a hydrocarbon in which one or more of the hydrogens has been replaced by nitro groups.
Heuristic:
  1. The molecule must contain at least one nitro group ([N+](=O)[O-]) attached directly to a carbon atom.
  2. Removing (i.e. “masking”) all atoms from nitro groups, the remaining heavy atoms (non-hydrogen) should be predominantly (>=90%)
     carbon atoms.
  3. For aromatic backbones: if more than one ring exists, the rings must be part of a fused system.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    
    A nitrohydrocarbon is defined as a hydrocarbon where one or more hydrogens has been replaced
    by nitro groups. In our heuristic we check that:
    
      1. The molecule contains at least one nitro group ([N+](=O)[O-]) attached directly to a carbon atom.
      2. When the nitro group atoms are removed, the remaining heavy-atom backbone is almost exclusively carbon (>=90%).
      3. If the backbone is aromatic and contains multiple rings, then the rings must be fused.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise.
        str: A reason describing the outcome.
    """
    # Parse the SMILES into an RDKit molecule.
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
        # In a nitro match, typically the first atom is the nitrogen.
        nitroN_idx = match[0]
        atomN = mol.GetAtomWithIdx(nitroN_idx)
        # Check neighbors of the nitro nitrogen.
        for neighbor in atomN.GetNeighbors():
            # The neighbor should be a carbon and not be part of the nitro group.
            if neighbor.GetAtomicNum() == 6 and (neighbor.GetIdx() not in match):
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
    
    # Create a modified copy of the molecule (the "backbone") by removing nitro group atoms.
    backbone = Chem.RWMol(mol)
    # Remove atoms in descending order to avoid index shift.
    for idx in sorted(nitro_atom_indices, reverse=True):
        try:
            backbone.RemoveAtom(idx)
        except Exception:
            pass
    backbone = backbone.GetMol()
    
    # Gather heavy atoms (exclude hydrogens) from the backbone.
    heavy_atoms = [atom for atom in backbone.GetAtoms() if atom.GetAtomicNum() != 1]
    if not heavy_atoms:
        return False, "Backbone is empty after removing nitro groups"
    
    # Count how many of these heavy atoms are carbon.
    carbon_count = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
    hydro_ratio = carbon_count / len(heavy_atoms)
    if hydro_ratio < 0.9:
        return False, f"Backbone has only {hydro_ratio*100:.1f}% carbons (requires at least 90%); extra heteroatoms present"
    
    # If the backbone contains aromatic rings, ensure that if multiple rings exist they form a fused system.
    ri = backbone.GetRingInfo()
    rings = ri.AtomRings()
    if len(rings) > 1:
        # For each ring, represent bonds as sets of frozensets (each representing a bond).
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
    
    return True, "Molecule is a nitrohydrocarbon: nitro group(s) are attached directly to a carbon and backbone is (mostly) hydrocarbon"

# Uncomment below for quick testing.
# test_smiles_list = [
#     "[O-][N+](=O)c1cccc2C3OC3C=Cc12",  # 1-Nitronaphthalene-5,6-oxide (true positive expected)
#     "[O-][N+](=O)c1ccc-2c(c1)-c1ccc([N+]([O-])=O)c3cccc-2c13",  # 3,9-Dinitrofluoranthene (true positive expected)
#     "CC(C)[N+]([O-])=O",  # 2-nitropropane (true positive expected)
#     "[O-][N+](=O)c1ccccc1",  # Nitrobenzene (true positive expected)
# ]
#
# for s in test_smiles_list:
#     result, reason = is_nitrohydrocarbon(s)
#     print(s, result, reason)