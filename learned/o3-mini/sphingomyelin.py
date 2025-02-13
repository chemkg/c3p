"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: Sphingomyelin
Definition:
  Any of a class of phospholipids in which the amino group of a sphingoid base is in amide linkage 
  with one of several fatty acids, while the terminal hydroxy group of the sphingoid base is 
  esterified to phosphorylcholine.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    A sphingomyelin has:
      - an amide bond linking a fatty acid to the sphingoid base (frequently indicated by a substructure such as C(=O)N[C@H] or C(=O)N[C@@H])
      - a phosphorylcholine head group esterified to the terminal hydroxyl of the sphingoid base (typically COP([O-])(=O)OCC[N+](C)(C)C)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as sphingomyelin, False otherwise.
        str: A reason describing the classification.
    """
    # Parse the SMILES to get the molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the phosphorylcholine head group.
    # This pattern targets: O-CH2-P(=O)(O-)OCH2CH2N+(C)(C)C.
    # The SMARTS below should match the common substructure seen in our examples.
    phosphocholine_smarts = "COP([O-])(=O)OCC[N+](C)(C)C"
    phospho_group = Chem.MolFromSmarts(phosphocholine_smarts)
    if not mol.HasSubstructMatch(phospho_group):
        return False, "Missing phosphorylcholine head group"
    
    # Check for a typical amide bond usually seen in sphingomyelin 
    # (the amino group of the sphingoid base linked to a fatty acid chain).
    # We check for both chiral representations.
    amide_smarts1 = "C(=O)N[C@H]"
    amide_smarts2 = "C(=O)N[C@@H]"
    amide1 = Chem.MolFromSmarts(amide_smarts1)
    amide2 = Chem.MolFromSmarts(amide_smarts2)
    if not (mol.HasSubstructMatch(amide1) or mol.HasSubstructMatch(amide2)):
        return False, "Missing amide linkage between the fatty acid and the sphingoid base"
    
    # Additional sanity checks:
    # 1. The molecule should have a sufficient number of carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Too few carbon atoms to be a sphingomyelin"
    
    # 2. Phospholipids like sphingomyelin generally have high molecular weights.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight is too low for sphingomyelin"
    
    return True, "Contains both a phosphocholine head group and an amide-linked fatty acid on a sphingoid base"

# Uncomment below lines to test the function with one of the provided SMILES examples:
# example_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCC(C)C"
# result, reason = is_sphingomyelin(example_smiles)
# print(result, reason)