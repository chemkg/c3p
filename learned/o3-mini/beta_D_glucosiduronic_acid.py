"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: beta-D-glucosiduronic acid conjugate
Definition: A glucosiduronic acid resulting from the formal condensation 
of any substance with beta-D-glucuronic acid to form a glycosidic bond.
The method identifies a pyranose glucuronic acid fragment (allowing for O‐ or N‐linkage)
and checks that the linking atom is attached to an external aglycone.
"""

from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) is a beta-D-glucosiduronic acid conjugate.
    The algorithm searches for a glucuronic acid moiety defined as a pyranose ring bearing a carboxyl group 
    and attached via a glycosidic bond through an oxygen or nitrogen. It then verifies that the linking atom 
    (the anomeric heteroatom) is connected to an external (non-sugar) fragment.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a beta-D-glucosiduronic acid conjugate, False otherwise.
        str: A message explaining the basis for the classification.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a relaxed SMARTS for the beta-D-glucuronic acid moiety.
    # We allow for the glycosidic linking atom to be either O or N:
    #   [O,N] - the linking heteroatom (i.e. the atom that bonds the sugar to the aglycone)
    #   C1 ... C(O1) - a six-membered (pyranose) ring in which one atom is the ring oxygen (specified by O1)
    #   and one of the ring carbons bears a carboxyl group (C(=O)O) appended to it.
    # This pattern is a relaxed approximation and does not enforce strict chiral configuration.
    sugar_smarts = "[O,N]C1C(O)C(O)C(O)C(O1)C(=O)O"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    if sugar_pattern is None:
        return None, "Error compiling SMARTS for glucuronic acid"
    
    # Search for all substructure matches (ignoring chirality for robustness).
    matches = mol.GetSubstructMatches(sugar_pattern, useChirality=False)
    if not matches:
        return False, "No beta-D-glucuronic acid substructure found"
        
    # For each identified sugar match, verify that the linking heteroatom provides a true glycosidic bond by 
    # having at least one neighbor that is not part of the sugar fragment.
    for match in matches:
        # match[0] corresponds to the atom matching the [O,N] (glycosidic linking atom).
        linking_atom = mol.GetAtomWithIdx(match[0])
        external_attachment = False
        for neighbor in linking_atom.GetNeighbors():
            # If the neighbor is not in the substructure match, we regard this as an external attachment.
            if neighbor.GetIdx() not in match:
                external_attachment = True
                break
        if external_attachment:
            return True, "Contains beta-D-glucuronic acid moiety attached via a glycosidic bond"
    
    # If we have found the sugar fragment but no external attachment (the linking atom only bonds within the sugar),
    # then the fragment is likely a free acid.
    return False, "Glucuronic acid moiety found but no glycosidic attachment to an external aglycone"