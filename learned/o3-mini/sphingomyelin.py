"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: Sphingomyelin

Definition:
  Sphingomyelins are phospholipids in which the amino group of a sphingoid base is in 
  amide linkage with a fatty acid and the terminal hydroxyl group of the sphingoid base is 
  esterified to phosphorylcholine. This implementation focuses on detecting a core motif:
  a fatty acid carbonyl (as a [CX3](=O)) bound directly to an amide nitrogen which is
  in turn immediately attached to a chiral carbon bearing a phosphorylcholine head group.
  
  We allow for slight variations in the phosphate representation (charged or neutral) as 
  well as in the chiral annotation.
"""

from rdkit import Chem

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    
    Our strategy is to first parse the molecule then to search for a specific structural
    “core” motif that captures a fatty acid carbonyl attached (as an amide) to a sphingoid 
    base that bears a phosphorylcholine head group. The SMARTS patterns are designed to require 
    that:
      - A carbonyl group ([CX3](=O)) is directly attached to a nitrogen,
      - That same nitrogen is directly attached to a chiral carbon (specified as [C@H] or [C@@H])
        that carries a phosphorylcholine moiety.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as sphingomyelin, False otherwise.
        str: A reason describing the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a set of SMARTS patterns that capture the sphingomyelin core.
    # The pattern looks for a carbonyl (fatty acid) directly attached to an amide nitrogen
    # which in turn is immediately connected to a chiral carbon substituted with a phosphorylcholine head.
    # We include variants for charged phosphate ([O-]) and neutral phosphate as well as for alternate chirality.
    core_patterns = [
        "[CX3](=O)N[C@H](COP(=O)(O)OCC[N+](C)(C)C)",    # neutral phosphate, chiral [C@H]
        "[CX3](=O)N[C@@H](COP(=O)(O)OCC[N+](C)(C)C)",   # neutral phosphate, chiral [C@@H]
        "[CX3](=O)N[C@H](COP([O-])(=O)OCC[N+](C)(C)C)",  # charged phosphate, chiral [C@H]
        "[CX3](=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)"  # charged phosphate, chiral [C@@H]
    ]
    
    # Try matching each core pattern.
    core_found = False
    for pattern in core_patterns:
        core_query = Chem.MolFromSmarts(pattern)
        if core_query is None:
            continue  # skip if SMARTS did not compile
        if mol.HasSubstructMatch(core_query):
            core_found = True
            break

    if not core_found:
        return False, ("Missing sphingomyelin core motif: expected a fatty acid carbonyl (C(=O)) "
                       "directly bound via an amide linkage to a chiral carbon carrying a "
                       "phosphorylcholine head group")
    
    # If the core motif is found, we classify the molecule as sphingomyelin.
    return True, ("Contains sphingomyelin core: fatty acid carbonyl amide directly linked to a "
                  "sphingoid base with phosphorylcholine head group")

# Example usage:
# test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C"
# result, reason = is_sphingomyelin(test_smiles)
# print(result, reason)