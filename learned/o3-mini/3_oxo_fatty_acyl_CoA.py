"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: 3-oxo-fatty acyl-CoA

An oxo fatty acyl-CoA results from the condensation of the thiol group of coenzyme A
with the carboxyl group of any 3-oxo fatty acid. A key motif is that the acyl chain
contains two carbonyl groups (one from the acid end and one as the 3-oxo substituent)
with a saturated (or branched) carbon in between, and it is attached via a thioester to
the CoA moiety.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    
    The function requires the molecule to have two characteristic fragments:
      1. A Coenzyme A moiety (detected by a fragment typical of the CoA “tail”).
      2. An acyl chain containing a thioester carbonyl linked to a 3-oxo group.
         For standard (open‐chain) cases, the minimal motif is the sequence
             –C(=O)–[CH2 or CH(R)]–C(=O)–S–
         We also allow one cyclic variant observed in some compounds.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple. First element is True if molecule seems to be a 
                     3-oxo-fatty acyl-CoA; otherwise False. The second element is 
                     a brief text description of the reason.
    """
    # Try to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the Coenzyme A moiety.
    # Many CoA derivatives contain the fragment "SCCNC(=O)CCNC(=O)" (which is part of the
    # pantetheine / adenosine portion). This is a heuristic and may not capture every case.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"

    # Look for the 3-oxo fatty acyl fragment.
    # For a standard 3-oxo acyl chain we expect a thioester with two carbonyl groups separated
    # by a saturated (sp3) carbon. The following SMARTS pattern covers cases like acetoacetyl:
    three_oxo_pattern = Chem.MolFromSmarts("C(=O)[CX4]C(=O)[S]")
    # In some examples the acyl chain is cyclic; for example, a cyclohexanone derivative may appear as:
    three_oxo_cyclic = Chem.MolFromSmarts("[S]C(=O)C1=CCCCC1=O")

    # Check if at least one of the patterns is found in the molecule
    has_three_oxo = mol.HasSubstructMatch(three_oxo_pattern) or mol.HasSubstructMatch(three_oxo_cyclic)
    if not has_three_oxo:
        return False, "3-oxo fatty acyl fragment not found"

    # Optional: one may add further checks on the length of the acyl chain or unsaturation
    # For example, many fatty acids should have a long alkyl chain; however, exceptions exist
    # (e.g. acetoacetyl-CoA or 2-methylacetoacetyl-CoA).
    
    # If all tests pass, classify as 3-oxo-fatty acyl-CoA
    return True, "Contains CoA moiety and 3-oxo fatty acyl fragment"

# Example usage (you can remove or comment these lines when using this as a module):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_3_oxo_fatty_acyl_CoA(test_smiles)
    print(result, reason)