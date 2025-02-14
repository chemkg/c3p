"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a 3-substituted propionyl-CoA(4-), False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for CoA(4-) core (phosphopantetheine part)
    coa_core_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC")
    if not mol.HasSubstructMatch(coa_core_pattern):
            return False, "CoA core substructure not found."

    # 2. Check for the thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
         return False, f"Found {len(thioester_matches)} thioester groups, need exactly 1"

    # 3. Check for propionyl group and the substituent at position 3, (and the connection to the thioester)
    #  This pattern matches a carbonyl carbon, with a connection to the sulfur, followed by 2 carbons with a substituent at the 3rd pos
    propionyl_pattern = Chem.MolFromSmarts("C(=O)C[C]")
    propionyl_matches = mol.GetSubstructMatches(propionyl_pattern)
    if len(propionyl_matches) != 1:
       return False, f"Found {len(propionyl_matches)} propionyl groups, need exactly 1"
    
    # Verify substitution on the beta carbon
    match = propionyl_matches[0]
    beta_carbon = mol.GetAtomWithIdx(match[2])
    if beta_carbon.GetDegree() < 3:
         return False, f"The beta carbon is not substituted"

    # 4. Verify the presence of the 4 negative charges, we look for 4 oxygen atoms that are deprotonated
    charged_o = 0
    for atom in mol.GetAtoms():
      if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
          charged_o += 1

    if charged_o != 4:
        return False, f"Incorrect number of negatively charged oxygens: {charged_o}, expected 4"


    return True, "Molecule matches the 3-substituted propionyl-CoA(4-) pattern."