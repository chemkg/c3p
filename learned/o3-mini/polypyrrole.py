"""
Classifies: CHEBI:38077 polypyrrole
"""
#!/usr/bin/env python
# polypyrrole.py
"""
Classifies: A compound composed of two or more pyrrole units.
A pyrrole unit is defined here as a five-membered aromatic ring 
that contains exactly one nitrogen and four carbon atoms.
"""

from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole, i.e. contains two or more pyrrole units,
    based on its SMILES string.
    
    A pyrrole ring is defined as a 5-membered ring in which:
      - All atoms are marked as aromatic.
      - There is exactly one nitrogen atom (atomic number 7).
      - The other four atoms are carbons (atomic number 6).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a polypyrrole, False otherwise.
        str: Detailed reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Note: Ensure that aromaticity is set (usually done automatically by RDKit).
    Chem.SanitizeMol(mol)
    
    # Obtain ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    pyrrole_unit_count = 0
    
    # Iterate over each ring in the molecule.
    for ring in atom_rings:
        # Check for 5-membered rings.
        if len(ring) != 5:
            continue
        
        # Get the atoms corresponding to the ring indices.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Check if every atom in the ring is aromatic.
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        
        # Count nitrogen atoms and verify that other atoms are carbon.
        n_nitrogen = 0
        non_carbon_present = False
        for atom in ring_atoms:
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 7:  # Nitrogen
                n_nitrogen += 1
            elif atomic_num != 6:  # Not carbon
                non_carbon_present = True
                break
        
        # For a pyrrole, we expect exactly 1 nitrogen and 4 carbons.
        if non_carbon_present:
            continue
        if n_nitrogen == 1:
            pyrrole_unit_count += 1
    
    # Classify according to the count of pyrrole units.
    if pyrrole_unit_count >= 2:
        return True, f"Contains {pyrrole_unit_count} pyrrole units (>=2); qualifies as a polypyrrole."
    else:
        return False, f"Found {pyrrole_unit_count} pyrrole unit(s); at least 2 are required for polypyrrole classification."

# Example usage (if running as a script):
if __name__ == "__main__":
    # You can test with a known polypyrrole-like molecule; for instance, porphyrins contain 4 pyrrole rings.
    test_smiles = "CC1=C(C)C2=C(C3=CC=CC=C3N2)N1"  # This is just a dummy SMILES for testing.
    result, reason = is_polypyrrole(test_smiles)
    print(result, reason)