"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: A disaccharide, defined as a compound in which two monosaccharides are joined by a glycosidic bond.
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is defined as two monosaccharide units (sugar rings) linked via a glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a disaccharide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Returns a tuple of tuples, each a collection of atom indices in a ring
    
    # Identify candidate sugar rings: rings of size 5 or 6 and containing exactly one oxygen.
    sugar_rings = []
    for ring in atom_rings:
        if len(ring) in (5, 6):  # pick rings of size typical for furanose or pyranose
            oxy_count = 0
            # Count oxygens in the ring
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    oxy_count += 1
            if oxy_count == 1:
                sugar_rings.append(set(ring))
    
    # For a disaccharide we expect exactly two sugar rings.
    if len(sugar_rings) != 2:
        return False, f"Expected 2 sugar rings, found {len(sugar_rings)} candidate(s)"
    
    # Check for a glycosidic bond: a bond connecting one atom in one sugar ring and one in the other.
    # Typically, glycosidic bonds involve an oxygen (O) atom bridging two monosaccharide units.
    ring1, ring2 = sugar_rings[0], sugar_rings[1]
    glyco_bond_found = False
    for bond in mol.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        # Check if the bond connects atoms from different sugar rings.
        in_ring1 = idx1 in ring1
        in_ring2 = idx1 in ring2
        in_ring1_2 = idx2 in ring1
        in_ring2_2 = idx2 in ring2
        if ((in_ring1 and idx2 in ring2) or (in_ring2 and idx2 in ring1)):
            # Further check if one of the atoms is an oxygen as typically expected for a glycosidic bond.
            atom1 = mol.GetAtomWithIdx(idx1)
            atom2 = mol.GetAtomWithIdx(idx2)
            if atom1.GetAtomicNum() == 8 or atom2.GetAtomicNum() == 8:
                glyco_bond_found = True
                break

    if not glyco_bond_found:
        return False, "No glycosidic bond connecting the two sugar rings found"
    
    return True, "Molecule contains two sugar rings linked by a glycosidic bond, consistent with a disaccharide"
    
# Example usage: Uncomment the following lines to test the function with a provided SMILES string.
# smiles_example = "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@H](O)[C@@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
# print(is_disaccharide(smiles_example))