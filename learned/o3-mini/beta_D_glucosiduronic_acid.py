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
    Determines if a molecule (given by its SMILES string) is a beta-D-glucosiduronic acid conjugate.
    The method works by searching for a beta-D-glucuronic acid moiety represented as a pyranose ring that
    carries an exocyclic carboxyl group and is attached via a glycosidic bond to an aglycone.
    
    To improve reliability the SMARTS pattern is written in a canonical form and chiral matching is disabled.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a beta-D-glucosiduronic acid conjugate, False otherwise.
        str: A reason message explaining the classification.
    """
    
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern representing beta-D-glucuronic acid.
    # This pattern shows the glycosidic linking oxygen (first atom) attached to a pyranose ring
    # that has a carboxyl (C(=O)O) substituent.
    # For matching purposes, we ignore chirality (useChirality=False).
    glucuronide_smarts = "O[C@@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1C(=O)O"
    sugar_pattern = Chem.MolFromSmarts(glucuronide_smarts)
    if sugar_pattern is None:
        return None, None  # SMARTS pattern failed to compile
    
    # Find matches for the glucuronic acid fragment while ignoring stereochemistry.
    matches = mol.GetSubstructMatches(sugar_pattern, useChirality=False)
    if not matches:
        return False, "No beta-D-glucuronic acid substructure found"
    
    # Verify that the linking oxygen (atom 0 in the pattern) actually forms a glycosidic bond
    # by checking if it is attached to an atom outside of the sugar fragment.
    for match in matches:
        glyco_O = mol.GetAtomWithIdx(match[0])
        external_attachment_found = False
        for neighbor in glyco_O.GetNeighbors():
            if neighbor.GetIdx() not in match:
                # To reduce false positives, we check that the external neighbor is not simply another oxygen.
                if neighbor.GetAtomicNum() != 8:
                    external_attachment_found = True
                    break
        if external_attachment_found:
            return True, "Contains beta-D-glucuronic acid moiety attached via a glycosidic bond"
            
    return False, "Glucuronic acid moiety found but no glycosidic attachment to an external aglycone"