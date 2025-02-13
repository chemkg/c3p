"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA(4-)
Definition:
  "An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups 
   of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3."

This classifier checks for:
  - A 3-hydroxy fatty acyl chain. We now require the presence of a hydroxylated stereocenter
    on a carbon that is directly attached to a CH2 group which in turn is attached to a thioester group.
    (SMARTS: "[C;H1]([OX2])[CH2]CC(=O)S" will match both (R) and (S) centers.)
  - A CoA moiety fragment as indicated by the substring "SCCNC(=O)CCNC(=O)".
  - At least 4 deprotonated oxygens in the phosphate/diphosphate parts. Since RDKit may not register
    a formal charge if not explicit, we also count occurrences of the literal substring "[O-]" in the SMILES.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    The molecule is required to contain:
      - A 3-hydroxy fatty acyl chain. In our classifier we look for a carbon with a hydroxyl 
        (either [C@H](O) or [C@@H](O)) immediately followed by a CH2 then CC(=O)S.
      - A CoA moiety fragment (we look for the fragment "SCCNC(=O)CCNC(=O)").
      - At least 4 deprotonated phosphate oxygen atoms. Because RDKit may not register 
        formal charges if not explicit, we also check the SMILES text for occurrences of "[O-]".
    
    Returns:
        (bool, str): True and an explanation if the molecule meets our criteria,
                     False and the reason it does not.
    """
    
    # Parse the input SMILES string to create an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for the 3-hydroxy fatty acyl chain.
    # Here we require a non-aromatic carbon bearing one hydrogen and an explicit O (hydroxyl),
    # followed by a CH2 and ultimately a thioester group.
    acyl_pattern = Chem.MolFromSmarts("[C;H1]([OX2])[CH2]CC(=O)S")
    if acyl_pattern is None:
        return False, "Failed to generate acyl substructure SMARTS pattern"
    
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No clear 3-hydroxy fatty acyl substructure ([C;H1](O)[CH2]CC(=O)S) found"

    # SMARTS for a CoA moiety fragment.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if coa_pattern is None:
        return False, "Failed to generate CoA moiety SMARTS pattern"
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety fragment (SCCNC(=O)CCNC(=O)) detected"

    # Count the number of oxygens with a formal charge of -1 from the molecule.
    neg_oxygens_rdkit = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    
    # In many of our provided examples the SMILES include explicit charges, but sometimes RDKit does not register them.
    # So we also count the literal occurrences of "[O-]" in the SMILES string.
    neg_oxygens_smiles = smiles.count("[O-]")
    
    # Use the maximum of the two counts.
    neg_oxygens = max(neg_oxygens_rdkit, neg_oxygens_smiles)
    
    if neg_oxygens < 4:
        return False, f"Insufficient deprotonated phosphate oxygens: found {neg_oxygens}, need at least 4"
    
    # All key features were found.
    return True, "Molecule contains the 3-hydroxy fatty acyl chain and CoA moiety with at least 4 deprotonated oxygens (CoA(4-))"

# Example usage (for testing, uncomment the following lines):
# test_smiles = "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC[C@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
# print(is_3_hydroxy_fatty_acyl_CoA_4__(test_smiles))