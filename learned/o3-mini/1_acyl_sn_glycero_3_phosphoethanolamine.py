"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
Definition: A 1-O-acylglycerophosphoethanolamine having (R)-configuration.
The molecule should have a fatty acyl chain esterified at the sn-1 position,
a phosphoethanolamine head group attached at the sn-3 position of glycerol,
and the chiral center (at the sn-1 carbon of the glycerol backbone) must be (R).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    The function checks for:
      1. A phosphoethanolamine head group.
      2. An acyl ester group attached to a glycerol backbone.
      3. A chiral center at the linkage that has (R)-configuration.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise.
        str: Reason for the classification decision.
    """

    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure that stereochemistry is assigned
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # 1. Check for the phosphoethanolamine head group.
    #   We look for a fragment like: –O–P(=O)(O)–O–C–C–N.
    #   In most examples the pattern "COP(O)(=O)OCCN" is present.
    head_smarts = "COP(O)(=O)OCCN"
    head_pattern = Chem.MolFromSmarts(head_smarts)
    if head_pattern is None:
        return False, "Invalid head group SMARTS pattern"
    if not mol.HasSubstructMatch(head_pattern):
        return False, "Phosphoethanolamine head group not found"

    # 2. Look for the 1-O-acyl ester group.
    # The acyl ester group linking an acyl chain to the glycerol backbone is expected to appear as:
    #   R-C(=O)O[C@H]... or R-C(=O)O[C@@H]...
    # We will try matching using two SMARTS patterns to cover both chiral notations.
    acylester_smarts1 = "C(=O)O[C@H]"
    acylester_smarts2 = "C(=O)O[C@@H]"
    acylester_pattern1 = Chem.MolFromSmarts(acylester_smarts1)
    acylester_pattern2 = Chem.MolFromSmarts(acylester_smarts2)
    matches1 = mol.GetSubstructMatches(acylester_pattern1) if acylester_pattern1 is not None else []
    matches2 = mol.GetSubstructMatches(acylester_pattern2) if acylester_pattern2 is not None else []
    acylester_matches = list(matches1) + list(matches2)
    if not acylester_matches:
        return False, "Acyl ester group (1-O-acyl) not found"

    # 3. Verify the stereochemistry of the glycerol chiral center.
    #    In our SMARTS the chiral center is the third atom in the match tuple (index 2).
    #    We then check its CIP code. We require it to be "R".
    valid_chiral = False
    for match in acylester_matches:
        # match is a tuple of indices; index 2 is the chiral carbon from our SMARTS.
        # Note: Depending on how the substructure is mapped, we check if this atom has a CIP code "R".
        chiral_idx = match[2]
        atom = mol.GetAtomWithIdx(chiral_idx)
        # Try to get the CIP code; if not set, then skip this match.
        try:
            cip = atom.GetProp('_CIPCode')
        except KeyError:
            cip = None
        if cip == "R":
            valid_chiral = True
            break
    if not valid_chiral:
        return False, "No acyl ester chiral center with (R)-configuration found"

    # Additional optional checks could be performed (e.g. ensuring the ester is linked to a long chain)
    # but for our current classification we assume these two criteria suffice.

    return True, "Molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine with (R)-configuration"

# Example usage:
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN"  # 1-stearoyl-sn-glycero-3-phosphoethanolamine
    result, reason = is_1_acyl_sn_glycero_3_phosphoethanolamine(test_smiles)
    print(result, reason)