"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: Anthoxanthin pigments – a subclass of flavonoids

Anthoxanthins are typically flavones or flavonols having a 2‐phenylchromen‐4‐one scaffold.
They have a fused system (rings A and C forming a benzopyran-4-one with a keto group at position 4)
with an attached benzene ring (ring B) at position 2. They are decorated with hydroxyl and/or methoxy groups
and are typically water soluble (i.e. low lipophilicity).

The classification criteria used here are:
  1. Valid SMILES string.
  2. Molecule must contain a 2‐phenylchromen‐4‐one (flavonoid) core.
     We identify this by a SMARTS pattern that captures the entire three‐ring scaffold.
     Our SMARTS "c1ccc(cc1)-c2oc3ccccc3c(=O)c2" is intended to match a substituted flavone core.
     (We further require that any core match covers exactly 15 atoms.)
  3. Presence of at least one free hydroxyl group ("[OX2H]").
  4. A low Crippen logP (<= 3.0) to support water‐solubility.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Crippen

def is_anthoxanthin(smiles: str):
    """
    Determines whether the molecule given as a SMILES string is likely to be an anthoxanthin.
    
    Criteria used:
      1. The molecule must be a valid structure.
      2. It must contain a 2‐phenylchromen‐4‐one core (flavonoid core).
         We search for a substructure matching the SMARTS:
           "c1ccc(cc1)-c2oc3ccccc3c(=O)c2"
         and further require that the found core consists of exactly 15 atoms.
      3. At least one free hydroxyl group ([OX2H]) must be present.
      4. The calculated Crippen logP must be <= 3.0.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where the first element is True if the molecule is classified as an anthoxanthin,
                     and False otherwise. The second element provides the reason.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the 2-phenylchromen-4-one flavonoid scaffold.
    # This pattern matches an aromatic benzene ring (ring B), attached to a benzopyran-4-one (rings A and C).
    # We expect a match covering 15 atoms in the core.
    flavonoid_core_smarts = "c1ccc(cc1)-c2oc3ccccc3c(=O)c2"
    flavonoid_core = Chem.MolFromSmarts(flavonoid_core_smarts)
    if flavonoid_core is None:
        return False, "Internal error: invalid SMARTS pattern"
    
    core_matches = mol.GetSubstructMatches(flavonoid_core)
    # Check that at least one match is found and that the matched set covers exactly 15 atoms.
    valid_core = False
    for match in core_matches:
        if len(match) == 15:
            valid_core = True
            break
    if not valid_core:
        return False, "No valid flavonoid core (2-phenylchromen-4-one) found"
    
    # Check for free hydroxyl groups (-OH); these are matched using "[OX2H]"
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    n_hydrox = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if n_hydrox < 1:
        return False, "Lacks free hydroxyl group(s) typical of anthoxanthins"
    
    # Calculate Crippen logP to assess water-solubility.
    logP = Crippen.MolLogP(mol)
    if logP > 3.0:
        return False, f"LogP too high ({logP:.2f}); molecule appears too lipophilic for a typical anthoxanthin"
    
    return True, "Molecule contains a valid 2-phenylchromen-4-one core, free hydroxyl group(s), and low logP, consistent with an anthoxanthin"

# Example usage:
# test_smiles = "COC1=CC=C(C=C1O)C1=CC(=O)C2=C(O)C=C(O)C=C2O1"  # diosmetin
# print(is_anthoxanthin(test_smiles))