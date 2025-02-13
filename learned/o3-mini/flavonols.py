"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: CHEBI: Flavonols
Definition: Any hydroxyflavone in which the ring hydrogen at position 3 of the heterocyclic (chromen) ring is replaced by a hydroxy group.
Our strategy:
  - Use a SMARTS pattern that represents the 3-hydroxyflavone core, i.e.
    2-phenyl-3-hydroxychromen-4-one.
  - The SMARTS below is based on the structure:
      O=c1c(O)cc2oc(-c3ccccc3)cc12
    which has a carbonyl (position 4), an –OH at the 3-position, the oxygen of the heterocycle, and a phenyl ring attached at position 2.
  - If the molecule contains that scaffold as a substructure (even if other substituents are attached), it is classified as a flavonol.
"""

from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a hydroxyflavone in which the ring hydrogen at position 3
    of the heterocyclic (chromen) ring is replaced by a hydroxy group (that is,
    the molecule contains a 3-hydroxyflavone scaffold).
    
    Our strategy is to directly search for a substructure corresponding to
    the 3-hydroxyflavone core using SMARTS.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a flavonol, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a 3-hydroxyflavone (flavonol) core.
    # The core structure is: O=c1c(O)cc2oc(-c3ccccc3)cc12
    #  • "O=c1" corresponds to a carbonyl group at position 4;
    #  • "c(O)" ensures that there is a hydroxyl group at position 3;
    #  • "oc(-c3ccccc3)cc12" defines the remainder of the fused ring system with the B-ring.
    flavonol_core_smarts = "O=c1c(O)cc2oc(-c3ccccc3)cc12"
    flavonol_pattern = Chem.MolFromSmarts(flavonol_core_smarts)
    if flavonol_pattern is None:
        return None, None  # something went wrong with SMARTS construction
    
    # Check if the molecule has the 3-hydroxyflavone substructure.
    if mol.HasSubstructMatch(flavonol_pattern):
        return True, "Molecule contains a flavonol (3-hydroxyflavone) scaffold."
    else:
        return False, "Molecule does not contain the required 3-hydroxyflavone scaffold for a flavonol."

# Example usage (you may remove or comment these out in production):
if __name__ == "__main__":
    # Some example SMILES strings (both from true positives and negatives)
    test_smiles = [
        # True positive: kaempferol 4'-O-beta-D-glucopyranoside
        "OC[C@H]1O[C@@H](Oc2ccc(cc2)-c2oc3cc(O)cc(O)c3c(=O)c2O)[C@H](O)[C@@H](O)[C@@H]1O",
        # True positive: 7,4'-dimethylkaempferol
        "C12=C(OC(C3=CC=C(OC)C=C3)=C(C1=O)O)C=C(OC)C=C2O",
        # Example flavone scaffold without 3-OH: (should not match as flavonol)
        "c1ccccc1-c2oc3ccccc3c(=O)c2"
    ]
    
    for smi in test_smiles:
        result, reason = is_flavonols(smi)
        print(f"SMILES: {smi}")
        print(f"Flavonol? {result} ({reason})\n")