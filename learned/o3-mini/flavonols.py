"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: Flavonols (CHEBI: Flavonols)
Definition: Any hydroxyflavone in which the ring hydrogen at position 3 of the heterocyclic (chromen) ring is replaced by a hydroxy group.
Our strategy:
  - Use a SMARTS pattern (based on 3-hydroxyflavone) corresponding to 2-phenyl-3-hydroxychromen-4-one.
  - To allow for derivatives (e.g. glycosides, extra substituents) and avoid failures because of stereochemistry,
    we remove stereochemical labels and set useChirality=False in the substructure search.
"""

from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol (i.e. contains a 3-hydroxyflavone scaffold)
    based on its SMILES string. A flavonol is defined as any hydroxyflavone in which the ring 
    hydrogen at position 3 of the heterocyclic (chromen) ring is replaced by a hydroxy group.
    
    Our strategy is to search for the 3-hydroxyflavone core using a SMARTS pattern.
    To improve matching (since prior attempts failed due to stereochemical and substituent issues)
    we remove stereochemistry and set useChirality=False.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a flavonol, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove stereochemical information to avoid mismatches due to chiral centers.
    Chem.RemoveStereochemistry(mol)
    
    # Define a SMARTS pattern for the 3-hydroxyflavone (flavonol) core.
    # The core structure is based on 2-phenyl-3-hydroxychromen-4-one.
    # SMARTS: O=c1c(O)cc2oc(-c3ccccc3)cc12
    # This pattern looks for:
    #   - a carbonyl (O=) at position 4,
    #   - a hydroxyl (OH) group at position 3,
    #   - and a connected aromatic (phenyl) B-ring at position 2.
    flavonol_core_smarts = "O=c1c(O)cc2oc(-c3ccccc3)cc12"
    flavonol_pattern = Chem.MolFromSmarts(flavonol_core_smarts)
    if flavonol_pattern is None:
        return None, None  # In case SMARTS construction fails
    
    # Perform the substructure search:
    # We use useChirality=False to ignore possible stereochemical differences.
    if mol.GetSubstructMatches(flavonol_pattern, useChirality=False):
        return True, "Molecule contains a flavonol (3-hydroxyflavone) scaffold."
    else:
        return False, "Molecule does not contain the required 3-hydroxyflavone scaffold for a flavonol."

# Example usage (you may comment these out for production):
if __name__ == "__main__":
    # Some test SMILES strings from the provided examples:
    test_smiles = [
        # kaempferol 4'-O-beta-D-glucopyranoside
        "OC[C@H]1O[C@@H](Oc2ccc(cc2)-c2oc3cc(O)cc(O)c3c(=O)c2O)[C@H](O)[C@@H](O)[C@@H]1O",
        # 7,4'-dimethylkaempferol
        "C12=C(OC(C3=CC=C(OC)C=C3)=C(C1=O)O)C=C(OC)C=C2O",
        # ombuin
        "COc1cc(O)c2c(c1)oc(-c1ccc(OC)c(O)c1)c(O)c2=O",
        # A flavone that lacks the 3-OH (should not be classified as flavonol)
        "c1ccccc1-c2oc3ccccc3c(=O)c2"
    ]
    
    for smi in test_smiles:
        result, reason = is_flavonols(smi)
        print(f"SMILES: {smi}")
        print(f"Flavonol? {result} ({reason})\n")