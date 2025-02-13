"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: Galactosylceramide – any cerebroside in which the monosaccharide head group is galactose.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is defined as a glycosphingolipid in which the ceramide (sphingoid base linked via an amide
    bond) bears a single galactose moiety.
    
    This function uses two substructure queries:
       1. A galactose sugar substructure (allowing for both α and β configurations).
       2. An amide bond (C(=O)N) in a molecular context that suggests a ceramide.
    Additionally, a molecular weight cutoff is employed to weed out too‐small molecules.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a galactosylceramide, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define SMARTS for the galactose ring.
    # Many galactosylceramides contain a hexopyranose with an exocyclic –CH2OH substituent.
    # Two variants are provided here to allow both alpha and beta linkages.
    beta_gal = Chem.MolFromSmiles("CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O")
    alpha_gal = Chem.MolFromSmiles("CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O")
    has_gal = False
    if beta_gal and mol.HasSubstructMatch(beta_gal):
        has_gal = True
    elif alpha_gal and mol.HasSubstructMatch(alpha_gal):
        has_gal = True
    if not has_gal:
        return False, "No galactose sugar head group found"

    # Check for the presence of an amide bond. The presence of "C(=O)N" with a long acyl fragment suggests a ceramide backbone.
    # (This is a simplified check.)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found – lacking a ceramide backbone"

    # Check overall molecular weight (a crude filter since typical sphingolipids are >500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low to be a glycosphingolipid or ceramide derivative"
    
    # Optionally, one could also check for the presence of a long aliphatic chain.
    # For instance, counting carbon atoms in the longest chain or looking for a minimum number of aliphatic carbons.
    # Here we count total carbons:
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbon atoms to represent a typical ceramide aliphatic chain"
    
    return True, "Contains a galactose sugar head group and an amide-linked ceramide backbone (galactosylceramide)"

# Example test (uncomment the following lines to test):
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_galactosylceramide(test_smiles)
# print(result, reason)