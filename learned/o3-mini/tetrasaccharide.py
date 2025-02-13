"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: Tetrasaccharide
Definition: An oligosaccharide comprising four monomeric monosaccharide units.
The classifier looks for exactly four sugar rings (typical 5- or 6-membered rings with one oxygen)
in a single connected molecule.
"""

from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determine if a molecule (given as a SMILES string) is a tetrasaccharide.
    It checks for:
      - Valid SMILES and one single (connected) component.
      - Exactly 4 sugar rings (either 5-membered (furanose) with 4 carbons +1 oxygen, 
        or 6-membered (pyranose) with 5 carbons +1 oxygen).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): True with reason if molecule qualifies as a tetrasaccharide;
                     otherwise, False with a reason.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Ensure the molecule is a single connected entity (tetrasaccharides are one molecule)
    mol_smiles = Chem.MolToSmiles(mol)
    if '.' in mol_smiles:
        return False, "Molecule contains multiple disconnected fragments"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Define a helper function to check if a ring looks like a carbohydrate ring.
    def is_sugar_ring(atom_indices):
        n_atoms = len(atom_indices)
        # For typical cyclic sugars: pyranose (6-membered) or furanose (5-membered)
        if n_atoms not in (5, 6):
            return False
        # Count the number of oxygens and carbons within the ring.
        oxygen_count = 0
        carbon_count = 0
        for idx in atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                carbon_count += 1
        # A typical pyranose should have 6 atoms: 1 oxygen and 5 carbons,
        # a furanose should have 5 atoms: 1 oxygen and 4 carbons.
        if n_atoms == 6 and oxygen_count == 1 and carbon_count == 5:
            return True
        if n_atoms == 5 and oxygen_count == 1 and carbon_count == 4:
            return True
        return False

    # Count the rings that qualify as carbohydrate rings.
    sugar_ring_count = 0
    for ring in ring_info:
        if is_sugar_ring(ring):
            sugar_ring_count += 1

    # Check if exactly four sugar rings were found.
    if sugar_ring_count != 4:
        return (False, f"Found {sugar_ring_count} candidate carbohydrate rings, need exactly 4")
    
    # (Optional additional checks could include verifying that the rings are connected via glycosidic bonds.)
    # For now, assume if the molecule is connected and has 4 sugar rings, it is a tetrasaccharide.
    
    return True, "Molecule contains exactly 4 monosaccharide (sugar) units"

# For testing you can run:
if __name__ == '__main__':
    test_smiles = "OC[C@H]1O[C@@H](OC[C@H]2O[C@@H](OC[C@H]3O[C@@H](OC[C@H]4OC(O)[C@H](O)[C@@H](O)[C@H]4O)[C@H](O)[C@@H](O)[C@H]3O)[C@H](O)[C@@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O"
    result, reason = is_tetrasaccharide(test_smiles)
    print(result, reason)