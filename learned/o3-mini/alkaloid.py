"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: Alkaloid
Definition: Any of the naturally occurring, basic nitrogen compounds (mostly heterocyclic) found
mostly in the plant kingdom (but also in bacteria, fungi, and animals). By extension, certain neutral
compounds biogenetically related to basic alkaloids are also classed as alkaloids. Compounds in which the
nitrogen is exocyclic (for example dopamine, mescaline, serotonin, etc.) are usually classed as amines rather than alkaloids.
This simple heuristic uses the presence of at least one nitrogen atom in a ring as an indication of an alkaloid.
"""
from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a given molecule (as a SMILES string) is likely an alkaloid.
    The heuristic used is:
      1) The SMILES must be parsable.
      2) The molecule must contain at least one nitrogen atom.
      3) At least one nitrogen should be embedded in a ring system (heterocyclic), because
         molecules with only exocyclic nitrogen functionality are usually considered simple amines.
      
    Note: This heuristic is very simplified. Real alkaloids are defined by their biogenetic origin
          and additional properties, and non-alkaloids such as amino acids, peptides, nucleotides, etc.
          are not usually considered as alkaloids.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule matches the alkaloid heuristic, False otherwise
        str: Explanation for the decision
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    total_nitrogen = 0
    ring_nitrogen = 0

    # Iterate over all atoms to count nitrogen and check ring membership
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atomic number is 7
            total_nitrogen += 1
            if atom.IsInRing():
                ring_nitrogen += 1

    # Must have at least one nitrogen atom
    if total_nitrogen == 0:
        return False, "No nitrogen atoms present, so unlikely to be an alkaloid"
    
    # Require at least one ring nitrogen; if all nitrogens are exocyclic, assume it's an amine
    if ring_nitrogen == 0:
        return False, "Nitrogen atoms present are exclusively exocyclic, suggesting an amine rather than an alkaloid"

    return True, f"Found {ring_nitrogen} ring nitrogen atom(s) out of {total_nitrogen} total nitrogen(s): consistent with an alkaloid classification"
    
# Example usage (can be removed or commented out)
if __name__ == "__main__":
    test_smiles = "CN1C2=C(C=CC(=C2)OC)C3=C1[C@@H](N(CC34CCN(CC4)C(=O)NC5=CC(=CC=C5)F)C(=O)CN6CCOCC6)CO"  # example from user list
    result, reason = is_alkaloid(test_smiles)
    print("Alkaloid:", result)
    print("Reason:", reason)