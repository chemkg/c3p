"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: beta-D-glucosiduronic acid conjugate
Definition: A glucosiduronic acid resulting from the formal condensation 
of any substance with beta-D-glucuronic acid to form a glycosidic bond.
The method identifies a pyranose glucuronide ring that carries a carboxyl group
and checks that the glycosidic oxygen (the first atom in the SMARTS) links 
to an external aglycone.
"""

from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) is a beta-D-glucosiduronic acid conjugate.
    The algorithm works by searching for a glucuronide ring fragment using a relaxed,
    non-chiral SMARTS pattern. This fragment is defined as a six-membered ring containing one oxygen
    and five carbons with one of the carbons substituted with a carboxyl group.
    Then we verify that the glycosidic oxygen (the first atom in the SMARTS match)
    is attached to an external aglycone.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a beta-D-glucosiduronic acid conjugate, False otherwise.
        str: A reason message explaining the classification.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a relaxed SMARTS pattern for the beta-D-glucuronide moiety.
    # Here, the pattern "OC1C(O)C(O)C(O)C(O1)C(=O)O" represents:
    #   - A glycosidic oxygen (first 'O') that links an external group to the sugar.
    #   - A six-membered (pyranose) ring: C1 ... O1 where the ring consists of one oxygen and five carbons.
    #   - A carboxyl group (C(=O)O) appended to one of the ring carbons.
    sugar_smarts = "OC1C(O)C(O)C(O)C(O1)C(=O)O"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    if sugar_pattern is None:
        return None, None  # SMARTS pattern failed to compile
    
    # Find substructure matches for the glucuronide ring.
    # We ignore chiral specifications for robustness.
    matches = mol.GetSubstructMatches(sugar_pattern, useChirality=False)
    if not matches:
        return False, "No beta-D-glucuronic acid substructure found"
    
    # For each identified sugar match, verify that the glycosidic oxygen (first atom in the pattern)
    # forms a glycosidic bond by having an external (non-sugar) connection.
    for match in matches:
        glyco_O = mol.GetAtomWithIdx(match[0])
        external_attachment = False
        # Check every neighbor of the glycosidic oxygen.
        for neighbor in glyco_O.GetNeighbors():
            # If the neighbor is not part of the matched sugar fragment, then it is considered an aglycone attachment.
            if neighbor.GetIdx() not in match:
                external_attachment = True
                break
        if external_attachment:
            return True, "Contains beta-D-glucuronic acid moiety attached via a glycosidic bond"
    
    # If we found the sugar fragment but the glycosidic oxygen is not linked outside,
    # then the fragment is likely a free acid.
    return False, "Glucuronic acid moiety found but no glycosidic attachment to an external aglycone"