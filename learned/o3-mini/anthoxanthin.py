"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: Anthoxanthin pigments (a type of flavonoid pigments, typically water‐soluble and ranging in color from white/colorless to yellow)
Examples include diosmetin, apigenin derivatives, luteolin derivatives, and many glycosylated forms.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Crippen

def is_anthoxanthin(smiles: str):
    """
    Determines if a given molecule (provided as a SMILES string) is likely to be an anthoxanthin.
    
    Anthoxanthins are a subclass of flavonoids (typically flavones or flavonols) that
    contain a flavonoid core – a 2-phenylchromen-4-one ring system – often decorated
    with hydroxyl (and/or methoxy) groups and sometimes glycosidic linkages. These pigments
    are water‐soluble (i.e. low lipophilicity) and their colors range from colorless/white to yellow.
    
    The heuristics used in this function include:
      1. Valid SMILES parsing.
      2. The presence of a flavonoid core (via a looser SMARTS pattern).
      3. The presence of at least one hydroxyl group.
      4. A minimum of 3 rings corresponding to rings A, B, and C.
      5. A relatively low Crippen MolLogP (water‐solubility).
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an anthoxanthin, False otherwise.
        str: Reason explaining the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a looser SMARTS pattern for the flavonoid core.
    # This pattern captures an aromatic benzene ring fused to a heterocyclic ring with a keto group.
    # It is based on the 2-phenylchromen-4-one scaffold, but allows for substituents.
    flavonoid_core = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)c(c2)")
    if not mol.HasSubstructMatch(flavonoid_core):
        return False, "No flavonoid core (2-phenylchromen-4-one) found"
    
    # Check for the presence of at least one hydroxyl group.
    hydroxyl = Chem.MolFromSmarts("[OX2H]")
    n_hydrox = len(mol.GetSubstructMatches(hydroxyl))
    if n_hydrox < 1:
        return False, "Lacks free hydroxyl groups typical of anthoxanthins"
    
    # Verify that the molecule has at least three rings (typical flavonoid structure has rings A, B, and C).
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Fewer than 3 rings; does not match the typical flavonoid architecture"
    
    # Calculate the Crippen logP value as a measure of lipophilicity.
    logP = Crippen.MolLogP(mol)
    # Anthoxanthins are water-soluble pigments. While the precise cutoff is uncertain,
    # we use a conservative threshold (for example, below 3.0) to allow some substitutions.
    if logP > 3.0:
        return False, f"LogP too high ({logP:.2f}); molecule appears too lipophilic for a typical anthoxanthin"
    
    return True, "Molecule contains a loose flavonoid core, adequate hydroxylation, three rings, and low logP, consistent with an anthoxanthin"