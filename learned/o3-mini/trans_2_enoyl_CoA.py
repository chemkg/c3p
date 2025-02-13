"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
An unsaturated fatty acyl-CoA that results from the formal condensation 
of the thiol group of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.
Improved version: attempts to restrict the classification to cases where the 
thioester group has an adjacent alpha atom that is engaged in a trans (E) double bond.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.

    For a positive classification, the molecule must have:
      1. A recognizable Coenzyme A moiety; here we look for the adenine fragment.
      2. At least one thioester group defined as [CX3](=O)[S;X2].
      3. In that thioester, the acyl (alpha) carbon (attached to the carbonyl group) must be
         double-bonded to a beta carbon with explicit trans (E) stereochemistry.
         (This implies the acyl chain is a 2-enoyl group.)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as trans-2-enoyl-CoA, False otherwise.
        str: Reason for classification outcome.
    """

    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure stereochemistry is properly assigned (to capture E/Z labels)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Check for presence of a Coenzyme A portion:
    # We search for an adenine/purine fragment found in Coenzyme A.
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A moiety (adenine/purine fragment not found)"

    # Look for a thioester group, using a SMARTS that requires a carbonyl [CX3](=O)
    # bonded to a divalent sulfur [S;X2]:
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[S;X2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester group required for acyl-CoA formation"

    # Check each thioester match to see if it has the proper trans-2-enoyl fragment.
    # In a proper 2-enoyl fragment, the carbonyl carbon (C(=O)) is attached (via a single bond)
    # to an alpha carbon that is further involved in a double bond to a beta carbon.
    found_trans_enoyl = False
    for match in thioester_matches:
        # match[0] is the carbonyl carbon and match[1] is the sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the alpha carbon: one of the neighbors of the carbonyl that is
        # not the carbonyl oxygen (which is double-bonded) nor the sulfur.
        alpha_atoms = []
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetIdx() == sulfur_idx:
                continue
            # Skip any oxygen (the carbonyl O will have atomic number 8)
            if neighbor.GetAtomicNum() == 8:
                continue
            alpha_atoms.append(neighbor)
        
        # We expect exactly one candidate for the acyl chain side.
        if len(alpha_atoms) != 1:
            continue
        alpha_atom = alpha_atoms[0]

        # Check that the alpha carbon is involved in a double bond (to a beta atom)
        # that is explicitly defined as having trans (E) stereochemistry.
        beta_double_found = False
        for bond in alpha_atom.GetBonds():
            # We look for a double bond
            if bond.GetBondType() != Chem.BondType.DOUBLE:
                continue
            # Ignore the double bond if it is to the carbonyl (should not occur in a thioester)
            other_atom = bond.GetOtherAtom(alpha_atom)
            if other_atom.GetIdx() == carbonyl_idx:
                continue
            # Check that this double bond is flagged as E (trans) stereochemistry.
            if bond.GetStereo() == Chem.BondStereo.STEREOE:
                beta_double_found = True
                break
        if beta_double_found:
            found_trans_enoyl = True
            break

    if not found_trans_enoyl:
        return False, "No trans (E) double bond at the alpha-beta position adjacent to the thioester group found"

    return True, "Contains thioester group with an alpha,beta-unsaturated (trans-2-enoyl) acyl chain and Coenzyme A moiety"

# Example usage (for testing; remove or comment out in production):
if __name__ == "__main__":
    test_smiles = "CCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_trans_2_enoyl_CoA(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)