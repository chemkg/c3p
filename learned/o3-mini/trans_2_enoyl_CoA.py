"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
An unsaturated fatty acyl-CoA that results from the formal condensation 
of the thiol group of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.

    This classifier checks for:
      1. The presence of a Coenzyme A portion using a characteristic SMARTS pattern.
      2. A thioester group (C(=O)S) as the junction between the acyl chain and CoA.
      3. An alpha,beta-unsaturated acyl chain with a double bond adjacent to the thioester
         that is specified with E (trans) stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as trans-2-enoyl-CoA, False otherwise.
        str: Reason for classification outcome.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure correct perception of stereochemistry from the SMILES input
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Check for the presence of the Coenzyme A moiety.
    # We use a SMARTS pattern that captures a characteristic fragment of CoA.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A moiety"
    
    # Look for the thioester group "C(=O)S" that is typical for an acyl-CoA.
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester group required for acyl-CoA formation"
    
    # For each thioester match, check if it is attached to an acyl chain bearing a trans (E) double bond.
    found_trans = False
    for match in thioester_matches:
        # match contains indices of atoms in the pattern: [carbonyl C, carbonyl O, S].
        carbonyl_idx = match[0]  # the carbon in C(=O)S
        atom_carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the "alpha" carbon attached to the carbonyl.
        # Exclude the oxygen and sulfur that define the thioester.
        alpha_candidates = []
        for neighbor in atom_carbonyl.GetNeighbors():
            if neighbor.GetIdx() not in match[1:]:
                alpha_candidates.append(neighbor)
        if not alpha_candidates:
            continue  # No alpha carbon found; try next match
        
        alpha = alpha_candidates[0]
        
        # Search for a double bond departing from the alpha carbon.
        # In a trans-2-enoyl group, the alpha carbon should be part of a C=C double bond.
        for bond in alpha.GetBonds():
            # Look for a double bond
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other_atom = bond.GetOtherAtom(alpha)
                # Skip if this double bond is back to the carbonyl carbon
                if other_atom.GetIdx() == atom_carbonyl.GetIdx():
                    continue
                # Check the stereochemistry of the double bond.
                # RDKit marks E stereochemistry as Chem.BondStereo.STEREOE.
                if bond.GetStereo() == Chem.BondStereo.STEREOE:
                    found_trans = True
                    break
        if found_trans:
            break

    if not found_trans:
        return False, "No trans (E) double bond adjacent to the thioester (C(=O)S) group found"
    
    return True, "Contains thioester connected to an alpha,beta-unsaturated (trans-2-enoyl) acyl chain with Coenzyme A"

# Example usage (you may remove or comment these lines when integrating into a larger project):
if __name__ == "__main__":
    test_smiles = "CCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_trans_2_enoyl_CoA(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)