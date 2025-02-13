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

    To be classified as trans-2-enoyl-CoA, the molecule must have:
      1. A recognizable Coenzyme A moiety (we look for the adenine/purine fragment).
      2. A thioester group ([CX3](=O)[SX1]) linking the acyl chain to CoA.
      3. An acyl chain with an alpha,beta-unsaturation at the thioester junction:
         the carbonyl carbon (C(=O)) must be bonded to an alpha-carbon that in turn is
         double-bonded (with trans stereochemistry) to a beta-carbon.
         
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as trans-2-enoyl-CoA, False otherwise.
        str: Reason for classification outcome.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reassign stereochemistry (to help read E/Z information)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Check for the presence of a Coenzyme A portion. 
    # We look for a purine fragment which is part of the adenosine unit.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A moiety (adenine/purine fragment not found)"
    
    # Find the thioester group; use a SMARTS for a carbonyl bound to a sulfur:
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX1]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester group required for acyl-CoA formation"
    
    # Look for an alpha,beta-unsaturated fragment at the thioester junction.
    # In a trans-2-enoyl unit: R-CH=CH-C(=O)-S-CoA, where the carbonyl is C1, 
    # the alpha carbon is C2, and the beta carbon is C3.
    found_trans = False
    for match in thioester_matches:
        # In the pattern "[CX3](=O)[SX1]", match[0] is the carbonyl carbon and match[1] is the sulfur.
        carbonyl_idx = match[0]
        atom_carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the alpha carbon: among the neighbors of the carbonyl, ignore the oxygen (double-bonded)
        # and the sulfur (part of the thioester). In a valid enoyl, there must be one non-O/S neighbor.
        alpha_candidates = []
        for neighbor in atom_carbonyl.GetNeighbors():
            atomic_num = neighbor.GetAtomicNum()
            if atomic_num in (8, 16):  # skip O and S
                continue
            alpha_candidates.append(neighbor)
        if len(alpha_candidates) != 1:
            # Either no clear alpha carbon or ambiguous situation; check next thioester
            continue
        alpha = alpha_candidates[0]
        
        # The alpha carbon should be connected by a double bond to a beta carbon.
        double_bond_found = False
        for bond in alpha.GetBonds():
            if bond.GetBondType() != Chem.BondType.DOUBLE:
                continue
            # Get the atom on the other end of this double bond
            beta = bond.GetOtherAtom(alpha)
            # If beta is the carbonyl carbon, skip (this is the C=O bond)
            if beta.GetIdx() == atom_carbonyl.GetIdx():
                continue
            # Check that the double bond carries stereochemical information and is trans (E).
            # RDKit marks E as BondStereo.STEREOE.
            if bond.GetStereo() == Chem.BondStereo.STEREOE:
                double_bond_found = True
                break
        if double_bond_found:
            found_trans = True
            break

    if not found_trans:
        return False, "No trans (E) double bond at the alpha-beta position adjacent to the thioester group found"
    
    return True, "Contains thioester connected to an alpha,beta-unsaturated (trans-2-enoyl) acyl chain with Coenzyme A"


# Example usage (remove or comment out for production integration):
if __name__ == "__main__":
    test_smiles = "CCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_trans_2_enoyl_CoA(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)