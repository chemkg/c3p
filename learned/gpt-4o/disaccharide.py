"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound with two monosaccharides joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a generic monosaccharide unit (pyranose/furanose)
    mono_pattern = Chem.MolFromSmarts('[C@H]1(O)[C@@H](O)[C@@H](O)[C@@H](O)[C@H](O)O1|[C@H]1(O)[C@@H](O)[C@@H](O)[C@H]O1')
    
    # Identify monosaccharide units in the molecule
    mono_matches = mol.GetSubstructMatches(mono_pattern)
    if len(mono_matches) < 2:
        return False, "Less than two monosaccharide units found"

    # Define a SMARTS pattern for the glycosidic bond (C-O-C linkage between two sugars using recursive SMARTS)
    glycosidic_bond_pattern = Chem.MolFromSmarts('[C@H]O[C@H]')
    
    # Identify glycosidic bonds in the molecule
    glyco_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glyco_matches) < 1:
        return False, "No glycosidic bonds found between sugar units"

    # Check if each glycosidic bond connects two different monosaccharides from the matches
    mono_indices = set()
    for match in mono_matches:
        mono_indices.update(match)
    
    linked_monosaccharides = 0
    for bond_match in glyco_matches:
        if bond_match[0] in mono_indices and bond_match[2] in mono_indices:
            linked_monosaccharides += 1

    # Verify exactly two sugar units with appropriate linkage
    if len(mono_matches) == 2 and linked_monosaccharides >= 1:
        return True, "Contains two monosaccharide units joined by a glycosidic bond"
    else:
        return False, "Incorrect number of monosaccharide units or insufficient linkages"

# Test examples
examples = [
    "O([C@H]1[C@H](O)[C@H](O)C(O[C@@H]1CO)O)[C@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO",
    "O([C@H]1[C@@H](O)[C@H](O[C@@H](O)[C@@H]1O)CO)[C@@H]2OC[C@H](O)[C@H](O)[C@H]2O"
]

for example in examples:
    result, reason = is_disaccharide(example)
    print(f"SMILES: {example} -> Is disaccharide? {result}. Reason: {reason}")