"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: Rotenoid – Members of the class of tetrahydrochromenochromene.
Definition: Rotenoids consist of a cis‐fused tetrahydrochromeno[3,4-b]chromene skeleton and its substituted derivatives.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Our approach uses a substructure SMARTS pattern that attempts to capture the essential 
    fused bicyclic (tetrahydrochromenochromene) skeleton common to rotenoids.
    
    Note: Due to the large substitution variability among rotenoids, this query is a heuristic.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule likely belongs to the rotenoid class, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the core rotenoid scaffold.
    # This pattern attempts to capture a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton.
    # The pattern: 
    #  O1[C@@H]2CC[C@H]1c1ccc(O)cc1OC2
    # – O1: an oxygen atom starting a ring;
    # – [C@@H]2: a chiral carbon initiating a second ring;
    # – CC[C@H]1: a short saturated chain that closes the first ring;
    # – c1ccc(O)cc1: an aromatic ring bearing a hydroxyl substituent;
    # – OC2: an oxygen that links back to the second ring closure.
    rotenoid_pattern = Chem.MolFromSmarts("O1[C@@H]2CC[C@H]1c1ccc(O)cc1OC2")
    if rotenoid_pattern is None:
        # In the unlikely event an error occurs in SMARTS parsing.
        return False, "Error in SMARTS pattern"
    
    # Search for the substructure match in the molecule.
    if mol.HasSubstructMatch(rotenoid_pattern):
        return True, "Molecule contains the tetracyclic fused tetrahydrochromeno[3,4-b]chromene skeleton characteristic of rotenoids"
    else:
        return False, "Molecule does not contain the rotenoid core substructure"

# If run as a script, you might include some test examples (commented out).
# Example:
# test_smiles = "O1C2C(O)(C=3C(OC2)=CC(OC)=C(OC)C3)C(=O)C4=C1C=C(OC)C=C4"  # 12α-Hydroxymunduserone from the example list
# result, reason = is_rotenoid(test_smiles)
# print(result, reason)