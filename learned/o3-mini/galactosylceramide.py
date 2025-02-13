"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: Galactosylceramide – any cerebroside in which the monosaccharide head group is galactose.
Improved version: In addition to checking for a galactose ring and an amide bond, this version 
verifies that at least one amide bond involves a long acyl fragment (using a SMARTS that requires 10+ CH2 in a row)
which is typical of ceramide backbones.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is defined as a glycosphingolipid where a ceramide (i.e. a sphingoid base
    linked via an amide bond to a long-chain fatty acid) bears a single galactose.
    
    This function uses several filters:
       1. It requires the molecule to contain a galactose ring (either alpha or beta).
       2. It requires the presence of an amide bond.
       3. It further requires that at least one amide bond is located on a fatty acyl fragment;
          that is, the carbonyl carbon is bonded to a chain of at least 10 consecutive CH2 groups.
       4. A crude overall molecular weight and carbon count filter is applied.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if the molecule is classified as a galactosylceramide, False otherwise.
       str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define substructure for a galactose ring.
    # Many galactosylceramides contain a hexopyranose with an exocyclic –CH2OH.
    # We provide two versions to allow both beta and alpha linkages.
    beta_gal = Chem.MolFromSmiles("CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O")
    alpha_gal = Chem.MolFromSmiles("CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O")
    has_gal = False
    if beta_gal and mol.HasSubstructMatch(beta_gal):
        has_gal = True
    elif alpha_gal and mol.HasSubstructMatch(alpha_gal):
        has_gal = True
    if not has_gal:
        return False, "No galactose sugar head group found"
    
    # Check for the presence of a generic amide bond.
    # (C(=O)N is very common, so we use a SMARTS that finds any carbonyl attached to a nitrogen.)
    generic_amide = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    if not mol.HasSubstructMatch(generic_amide):
        return False, "No amide bond found – lacking a ceramide backbone"
    
    # Now check for an amide in which the acyl side is long (i.e. a fatty acyl chain).
    # We use a SMARTS pattern that requires a carbonyl group (C(=O)) followed by at least ten CH2 groups.
    # This pattern is a heuristic intended to detect an extended, unbranched alkyl chain.
    # Note: SMARTS quantifiers are supported in RDKit.
    fatty_acyl_smarts = "C(=O)[CH2]{10,}"
    fatty_acyl_pattern = Chem.MolFromSmarts(fatty_acyl_smarts)
    if not fatty_acyl_pattern:
        # Fallback in case SMARTS did not compile.
        return False, "Internal error: fatty acyl SMARTS pattern could not be parsed"
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No long fatty acyl chain detected on the amide – not a typical ceramide"
    
    # Additional crude filter: require molecular weight greater than 500 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low to be a glycosphingolipid or ceramide derivative"
    
    # Also ensure that the total carbon count is sufficiently high (at least 20 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbon atoms to represent a typical ceramide aliphatic chain"
    
    return True, ("Contains a galactose sugar head group and an amide-linked ceramide backbone "
                  "with a long fatty acyl chain (galactosylceramide)")

# Example test (uncomment to run):
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_galactosylceramide(test_smiles)
# print(result, reason)

# You may try additional test SMILES strings from the provided examples.