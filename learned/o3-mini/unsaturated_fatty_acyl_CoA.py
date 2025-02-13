"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: Unsaturated fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any unsaturated fatty acid.
This function checks for a thioester group (indicating linkage between acyl and CoA), a CoA moiety (using an adenine SMARTS pattern),
and finally the presence of at least one nonaromatic carbon–carbon double bond (indicating unsaturation in the fatty acid part).
"""

from rdkit import Chem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    
    The criteria used are:
      1. The molecule must be parsed successfully.
      2. It must contain a thioester group, i.e. a carbonyl (C=O) attached to a sulfur (S), 
         which is the signature of the linkage to CoA.
      3. It must contain a CoA moiety.
         Here we look for a characteristic purine/adenine ring pattern found in CoA.
      4. The fatty acyl chain must be unsaturated.
         We count the number of nonaromatic C=C bonds in the molecule.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is an unsaturated fatty acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for thioester group: look for a 
    # carbon atom double-bonded to an oxygen and singly bonded to a sulfur.
    #
    # SMARTS "[C](=O)[S]" will match the C(=O)S group characteristic of a thioester.
    thioester_pattern = Chem.MolFromSmarts("[C](=O)[S]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester group linking the fatty acyl part to CoA"
    
    # 2. Check for CoA moiety using a characteristic adenine/purine substructure.
    # In many CoA SMILES the fragment "n1cnc2c(ncnc12)" is present.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA moiety"
    
    # 3. Count nonaromatic C=C bonds which signal unsaturation.
    unsat_bond_count = 0
    for bond in mol.GetBonds():
        # Check if bond is a double bond, not aromatic and between carbon atoms.
        if (bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
            not bond.GetIsAromatic() and 
            bond.GetBeginAtom().GetAtomicNum() == 6 and 
            bond.GetEndAtom().GetAtomicNum() == 6):
            unsat_bond_count += 1

    if unsat_bond_count < 1:
        return False, "Fatty acyl chain appears to be saturated (no C=C bonds detected)"
    
    return True, ("Found a thioester linkage to CoA and at least one nonaromatic C=C double bond in the acyl chain, "
                  "consistent with an unsaturated fatty acyl-CoA")