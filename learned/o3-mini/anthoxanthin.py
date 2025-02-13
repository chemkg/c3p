"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: Anthoxanthin pigments (a type of flavonoid pigments, typically water-soluble and ranging in color from white/colorless to yellow)
Examples include diosmetin, apigenin derivatives, luteolin derivatives, and many glycosylated forms.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Crippen

def is_anthoxanthin(smiles: str):
    """
    Determines if a given molecule (via its SMILES string) is likely to be an anthoxanthin.
    
    Anthoxanthins are a subclass of flavonoids (often flavones or flavonols) with a characteristic
    2-phenylchromen-4-one (flavonoid) skeleton, generally decorated with hydroxyl groups and (often) glycosides,
    and are water-soluble pigments (white/colorless to yellow).
    
    The function uses several heuristic checks:
    1. The SMILES is parsed to a molecule.
    2. The molecule must contain the 2-phenylchromen-4-one flavonoid core.
    3. The molecule should have at least one free hydroxyl group.
    4. The molecule should contain at least three rings.
    5. The calculated logP (a measure of lipophilicity) is expected to be low for water‐soluble pigments.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an anthoxanthin, False otherwise.
        str: A message explaining the reasoning.
    """
    
    # Try to generate an RDKit molecule from the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the flavonoid core (2-phenylchromen-4-one scaffold).
    # This pattern is defined as a benzene ring attached at position 2 of a chromen-4-one.
    flavonoid_core = Chem.MolFromSmarts("c1ccccc1-c2c(=O)oc3ccccc23")
    if not mol.HasSubstructMatch(flavonoid_core):
        return False, "No flavonoid core (2-phenylchromen-4-one) found"
    
    # Check for free hydroxyl groups.
    # The [OX2H] SMARTS will match oxygen atoms with two connections (sp3) bonded to a hydrogen.
    hydroxyl = Chem.MolFromSmarts("[OX2H]")
    n_hydrox = len(mol.GetSubstructMatches(hydroxyl))
    if n_hydrox < 1:
        return False, "Lacks hydroxyl groups typical of anthoxanthins"
    
    # Verify that the molecule has three rings (the classical flavonoid structure has rings A, B, and C)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Fewer than 3 rings; does not match typical flavonoid architecture"
    
    # Calculate lipophilicity using Crippen's MolLogP.
    logP = Crippen.MolLogP(mol)
    # Anthoxanthins are water-soluble (generally low logP). Adjust the cutoff if needed.
    if logP > 2.5:
        return False, f"LogP too high ({logP:.2f}); molecule appears too lipophilic for typical anthoxanthins"
    
    return True, "Molecule contains a flavonoid core, adequate hydroxylation, three rings, and low logP, consistent with an anthoxanthin"