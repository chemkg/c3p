"""
Classifies: CHEBI:38077 polypyrrole
"""
#!/usr/bin/env python
"""
Classifies: A compound composed of two or more pyrrole units.
The definition here is that a pyrrole unit is a five‐membered aromatic ring 
that contains exactly one nitrogen atom (which may be substituted) and four carbon atoms.
In addition to finding at least two pyrrole units, we require that a significant fraction
of the molecule’s heavy atoms (i.e. non‐hydrogens) are in those pyrrole rings.
This extra condition aims to reduce mis‐classification of compounds with isolated pyrrole-like substructures.
"""

from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole; i.e. it is composed largely of pyrrole units,
    based on its SMILES string.
    
    Here a pyrrole unit is defined as a five-membered aromatic ring in which:
      - The ring is aromatic.
      - Exactly one atom is a nitrogen (atomic number 7), and the remaining four atoms are carbons (atomic number 6).
    
    In order to avoid false positives (e.g. compounds that merely contain a small pyrrole fragment),
    we also measure what fraction of the heavy (non‐hydrogen) atoms are part of the pyrrole matches.
    Only if the molecule has at least 2 pyrrole matches and if the union of all atoms in these matches 
    represents a sufficiently high fraction (here, 40%) of the molecule’s heavy atoms, 
    will it be classified as a polypyrrole.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a polypyrrole, False otherwise.
        str: Detailed reason for the classification decision.
    """
    
    # Attempt to parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize (this sets aromaticity, etc.).
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Count heavy atoms.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not heavy_atoms:
        return False, "No heavy atoms found in molecule"
    n_heavy = len(heavy_atoms)
    
    # Define a SMARTS that matches a pyrrole unit:
    # [n;a,r] means an aromatic nitrogen in a ring;
    # followed by four aromatic carbon atoms [c;a,r]. 
    # The ring closure "1" ensures a five-membered ring.
    pyrrole_smarts = "[n;a,r]1[c;a,r][c;a,r][c;a,r][c;a,r]1"
    pyrrole_query = Chem.MolFromSmarts(pyrrole_smarts)
    if pyrrole_query is None:
        return None, None  # Cannot compile SMARTS; unlikely
    
    # Find all pyrrole substructure matches.
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_query, uniquify=True)
    match_count = len(pyrrole_matches)
    
    # If less than 2 pyrrole matches, then we do not have a polypyrrole.
    if match_count < 2:
        return False, f"Found {match_count} pyrrole unit(s); at least 2 are required for polypyrrole classification."
    
    # Now, accumulate all atom indices (even if overlapping) that occur in any pyrrole match.
    pyrrole_atoms = set()
    for match in pyrrole_matches:
        pyrrole_atoms.update(match)
    
    # Compute the fraction of heavy atoms that are part of a pyrrole ring:
    fraction = len(pyrrole_atoms) / n_heavy
    
    # We choose a threshold that requires a significant portion of the molecule to be composed of pyrrole atoms.
    threshold = 0.4  # 40% of heavy atoms must be in matched pyrrole rings
    if fraction < threshold:
        return False, (f"Found {match_count} pyrrole unit(s), but only {fraction:.2f} "
                       f"of heavy atoms are in pyrrole rings (threshold {threshold}); "
                       "molecule is not predominantly composed of pyrrole units.")
    
    # If both conditions are satisfied, we classify as polypyrrole.
    return True, (f"Contains {match_count} pyrrole unit(s) and {fraction:.2f} fraction of heavy atoms are in pyrrole rings; "
                  "qualifies as a polypyrrole.")

# Example usage:
if __name__ == "__main__":
    # Test with a known polypyrrole-like molecule (e.g. a porphyrin which has 4 pyrrole units).
    test_smiles = "C1(=CC=2NC(=C(C2C)CCC(=O)O)C=C3NC(C(=C3CCC(=O)O)C)=O)NC(C(=C1C)C=C)=O"  # biotripyrrin-a (example)
    result, reason = is_polypyrrole(test_smiles)
    print(result, reason)