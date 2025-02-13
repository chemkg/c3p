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
      1. A recognizable Coenzyme A moiety: we look for the adenine fragment of adenosine.
      2. A thioester group, i.e. a carbonyl connected to a sulfur atom.
         (Note: we require the sulfur to be divalent: [CX3](=O)[S;X2])
      3. An acyl chain that is alpha,beta-unsaturated at the thioester junction:
         that is, the carbonyl carbon’s alpha carbon should be double-bonded 
         (with trans stereochemistry) to a beta carbon.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as trans-2-enoyl-CoA, False otherwise.
        str: Reason for classification outcome.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Optionally assign stereochemistry (to capture E/Z labels)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Check for the presence of a Coenzyme A portion (look for the adenine fragment)
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A moiety (adenine/purine fragment not found)"
    
    # Find the thioester group using a SMARTS pattern.
    # IMPORTANT: We use "[CX3](=O)[S;X2]" instead of "[CX3](=O)[SX1]" because the sulfur in a thioester is divalent.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[S;X2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester group required for acyl-CoA formation"
    
    # For each found thioester group, check the acyl chain unsaturation
    found_trans = False
    for match in thioester_matches:
        # match[0] is the carbonyl carbon in the thioester, and match[1] is the sulfur
        carbonyl_idx = match[0]
        atom_carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the alpha carbon. It must be one of the neighbors of the carbonyl
        # (excluding the oxygen of the carbonyl and the sulfur).
        alpha_candidates = []
        for neighbor in atom_carbonyl.GetNeighbors():
            atomic_num = neighbor.GetAtomicNum()
            if atomic_num in (8, 16):  # Skip oxygen (8) and sulfur (16)
                continue
            alpha_candidates.append(neighbor)
        
        # We expect exactly one candidate for the alpha carbon.
        if len(alpha_candidates) != 1:
            continue
        alpha_atom = alpha_candidates[0]
        
        # The alpha carbon should be connected by a double bond to a beta carbon.
        double_bond_found = False
        for bond in alpha_atom.GetBonds():
            if bond.GetBondType() != Chem.BondType.DOUBLE:
                continue
            # Identify the atom on the other side of the bond.
            beta_atom = bond.GetOtherAtom(alpha_atom)
            # Skip the bond back to the carbonyl carbon if it is present.
            if beta_atom.GetIdx() == atom_carbonyl.GetIdx():
                continue
            # Check that the double bond is marked as trans (E).
            if bond.GetStereo() == Chem.BondStereo.STEREOE:
                double_bond_found = True
                break
        if double_bond_found:
            found_trans = True
            break

    if not found_trans:
        return False, "No trans (E) double bond at the alpha-beta position adjacent to the thioester group found"
    
    return True, "Contains thioester group with an alpha,beta-unsaturated (trans-2-enoyl) acyl chain and Coenzyme A moiety"

# Example usage (remove or comment out when integrating in production):
if __name__ == "__main__":
    test_smiles = "CCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_trans_2_enoyl_CoA(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)