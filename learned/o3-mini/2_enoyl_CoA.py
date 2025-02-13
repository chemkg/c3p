"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: 2-enoyl-CoA
Definition: An unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if the molecule defined by the SMILES string is a 2-enoyl-CoA.
    The definition requires that the molecule has:
      1. A CoA moiety. For our purposes we check for a characteristic substructure of CoA.
      2. A thioester functional group where the acyl chain is unsaturated with a double bond between its C2 and C3.
         In an acyl chain, the carbonyl carbon is position 1 so the double bond must occur between the next two carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule qualifies as a 2-enoyl-CoA, otherwise False.
        str: Explanation of the classification result.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the unsaturated thioester group.
    # The pattern "[#6]=[#6]C(=O)S" looks for a carbon–carbon double bond 
    # followed by a carbonyl (C(=O)) that is directly attached to a sulfur.
    # This matches an acyl chain with a double bond between positions 2 and 3.
    enoyl_pattern = Chem.MolFromSmarts("[#6]=[#6]C(=O)S")
    if enoyl_pattern is None:
        return False, "Error in SMARTS for enoyl pattern"
    
    # Check if the molecule contains the alpha,beta unsaturated thioester fragment.
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No unsaturated thioester group with a double bond between positions 2 and 3 was detected."
    
    # Define a SMARTS pattern for a portion of the CoA moiety.
    # Many acyl-CoA structures share a common fragment such as "SCCNC(=O)CCNC(=O)"
    # which corresponds to part of the pantetheine and adenosine diphosphate units.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if coa_pattern is None:
        return False, "Error in SMARTS for CoA pattern"
    
    # Check if the molecule contains the CoA-related fragment.
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety fragment not detected."
    
    # Optionally, you can also check molecular descriptors (e.g., molecular weight) 
    # or other structural features typical for acyl-CoAs, but for the purpose 
    # of this classification our substructure matches are sufficient.
    
    return True, "Molecule contains a 2-enoyl-CoA fragment (unsaturated acyl thioester with CoA moiety)."

# Example usage (uncomment to test):
# test_smiles = "CCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
# result, reason = is_2_enoyl_CoA(test_smiles)
# print(result, reason)