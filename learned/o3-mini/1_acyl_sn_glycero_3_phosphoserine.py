"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoserine
Definition: An sn-glycerophosphoserine compound having an acyl substituent at the 1-hydroxy position.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    The molecule must contain both:
      - A phosphoserine headgroup (represented by a serine pattern OC[C](N)C(=O)O)
      - An acylated glycerol backbone at the sn-1 position (having an ester linkage O-C(=O)R)
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule fits the class, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the serine headgroup.
    # This pattern covers the serine portion generally depicted as OC[C@H](N)C(O)=O (ignoring strict chirality)
    serine_pattern = Chem.MolFromSmarts("OCC(N)C(=O)O")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "Serine headgroup (OC[C](N)C(=O)O) was not found"

    # Define a SMARTS pattern for the acylated glycerol moiety.
    # The pattern "OCC(O)COC(=O)[#6]" represents:
    #   - A glycerol backbone fragment: OCC(O)C...
    #   - An ester linkage: ...OC(=O)[#6]
    # The requirement of a carbon ([#6]) following the carbonyl ensures that an acyl chain is present.
    glycerol_acyl_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)[#6]")
    if not mol.HasSubstructMatch(glycerol_acyl_pattern):
        return False, "Acylated glycerol backbone (1-acyl substitution at the 1-hydroxy position) was not found"

    # If both patterns are present, then the molecule is classified in this chemical entity class.
    return True, ("Contains the serine headgroup and an acylated glycerol fragment, "
                  "consistent with 1-acyl-sn-glycero-3-phosphoserine")