"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: Sphingomyelin

Definition:
  Sphingomyelin is defined as any phospholipid in which the amino group of a sphingoid base is 
  in amide linkage with a fatty acid and the terminal hydroxyl group of the sphingoid base is 
  esterified to a phosphorylcholine. This program looks for a sphingomyelin “core” motif that 
  binds these features together.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    In this implementation we demand that the molecule:
      - Contains exactly one phosphorus atom and exactly two nitrogen atoms,
      - Has a sufficiently long carbon chain (>=30 carbons) and a molecular weight of at least 500 Da,
      - Has a core structural motif that links a fatty acid carbonyl (C(=O)) with an amide nitrogen 
        that is directly attached (with chiral specification) to a carbon bearing a phosphorylcholine head group.
    We allow for slight variations in the phosphate representation (charged or neutral).

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as sphingomyelin, False otherwise.
        str: A reason describing the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanity checks: count atoms of key elements.
    numP = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    numN = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    numC = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if numP != 1:
        return False, f"Expected exactly 1 phosphorus atom, but found {numP}"
    if numN != 2:
        return False, f"Expected exactly 2 nitrogen atoms, but found {numN}"
    if numC < 30:
        return False, f"Too few carbon atoms ({numC}) to be a sphingomyelin"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for sphingomyelin"
    
    # To better capture the sphingomyelin core, we look for a substructure that includes the fatty acid carbonyl
    # (as [CX3](=O)) connected to an amide nitrogen, which in turn is directly attached (with chiral specification)
    # to a carbon bearing the phosphorylcholine head group.
    # We allow variations in the phosphate pattern: sometimes it is drawn with an explicit [O-] and sometimes not.
    core_patterns = [
        "[CX3](=O)N[C@H](COP([O-])(=O)OCC[N+](C)(C)C)",  # charged phosphate, chiral as [C@H]
        "[CX3](=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)", # charged phosphate, alternate chirality
        "[CX3](=O)N[C@H](COP(=O)(O)OCC[N+](C)(C)C)",       # neutral phosphate, chiral as [C@H]
        "[CX3](=O)N[C@@H](COP(=O)(O)OCC[N+](C)(C)C)"       # neutral phosphate, alternate chirality
    ]
    core_found = False
    for pattern in core_patterns:
        core_mol = Chem.MolFromSmarts(pattern)
        if core_mol is None:
            continue
        if mol.HasSubstructMatch(core_mol):
            core_found = True
            break
    if not core_found:
        return False, "Missing sphingomyelin core structural motif (amide-linked fatty acid and phosphorylcholine on a sphingoid base)"
    
    return True, "Contains sphingomyelin core: fatty acid carbonyl amide linked to a sphingoid base with phosphorylcholine head group"

# Example (uncomment to test):
# example_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C"
# result, reason = is_sphingomyelin(example_smiles)
# print(result, reason)