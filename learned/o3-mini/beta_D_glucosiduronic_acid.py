"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: beta-D-glucosiduronic acid conjugate
Definition: A glucosiduronic acid resulting from the formal condensation 
of any substance with beta-D-glucuronic acid to form a glycosidic bond.
"""

from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule (given by SMILES string) belongs to the class 
    beta-D-glucosiduronic acid conjugates.
    
    This is defined by the presence of a beta-D-glucuronic acid moiety attached 
    via a glycosidic bond to an aglycone. The sugar unit is recognized by a 
    SMARTS pattern that matches a six-membered ring with a terminal carboxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a beta-D-glucosiduronic acid conjugate, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the beta-D-glucuronic acid moiety.
    # The pattern has a glycosidic oxygen (first atom) linking to the sugar ring.
    # It is assumed that the sugar has:
    #   - a pyranose ring with proper chiral descriptors as often seen for beta-D-glucuronic acid,
    #   - a carboxylic acid group (C(=O)O) attached at the ring.
    beta_glucuronide_smarts = "O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1C(=O)O)"
    sugar_pattern = Chem.MolFromSmarts(beta_glucuronide_smarts)
    if sugar_pattern is None:
        return None, None  # In case the SMARTS failed to compile
    
    # Find substructure matches for the glucuronic acid fragment
    matches = mol.GetSubstructMatches(sugar_pattern)
    if not matches:
        return False, "No beta-D-glucuronic acid substructure found"
    
    # Check if the glycosidic oxygen (the first atom in the SMARTS pattern) in at least one match
    # is linked to an atom outside of the matched fragment.
    for match in matches:
        # match[0] corresponds to the linking oxygen in our SMARTS.
        glyco_O = mol.GetAtomWithIdx(match[0])
        external_attachment_found = False
        # Iterate over neighbors of the glycosidic oxygen
        for neighbor in glyco_O.GetNeighbors():
            # If the neighbor atom is not part of the matched sugar fragment, it is from an aglycone.
            if neighbor.GetIdx() not in match:
                external_attachment_found = True
                break
        if external_attachment_found:
            return True, "Contains beta-D-glucuronic acid moiety attached via a glycosidic bond"
    
    return False, "Glucuronic acid moiety found but no glycosidic attachment to an external group"