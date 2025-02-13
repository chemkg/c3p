"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: 3-oxo-fatty acyl-CoA

An oxo fatty acyl-CoA results from the condensation of the thiol group of coenzyme A
with the carboxyl group of any 3-oxo fatty acid. A key motif is that the acyl chain
contains two carbonyl groups (one at the acid end and one as the 3-oxo substituent)
separated by one saturated (or branched) carbon, and it is attached via a thioester bond
to the CoA moiety.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    
    The method performs three main steps:
      1. Check that the SMILES parses and the overall formal charge is zero.
      2. Look for a CoA fragment (heuristically using the fragment "SCCNC(=O)CCNC(=O)").
      3. Look for a 3-oxo acyl fragment in which a thioester bond connects a 
         –C(=O)[CX4]C(=O)[S]– motif. In addition, we verify that the sulfur atom in
         the acyl fragment is directly connected to an atom that is part of the CoA fragment.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        (bool, str): Tuple. First element is True if the molecule is classified as a
                     3-oxo-fatty acyl-CoA; otherwise False. The second element is a
                     brief description of the reason.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule is neutral.
    if mol.GetFormalCharge() != 0:
        return False, "Molecule carries a non-zero formal charge"
    
    # Look for a Coenzyme A moiety.
    # Here we use a fragment common to many CoA derivatives.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "Coenzyme A moiety not found"
    # Create a set of all atom indices that are part of any CoA match
    coa_atoms = set()
    for match in coa_matches:
        for idx in match:
            coa_atoms.add(idx)
    
    # Define the standard 3-oxo pattern for an open-chain acyl fragment:
    # This minimal motif covers –C(=O)–[CH2 or CH(R)]–C(=O)–S–
    three_oxo_pattern = Chem.MolFromSmarts("C(=O)[CX4]C(=O)[S]")
    # And allow also one cyclic variant (e.g. a cyclohexanone derivative attached via thioester)
    three_oxo_cyclic = Chem.MolFromSmarts("[S]C(=O)C1=CCCCC1=O")
    
    # Check for the three_oxo matches (either standard or cyclic)
    matches_standard = mol.GetSubstructMatches(three_oxo_pattern)
    matches_cyclic = mol.GetSubstructMatches(three_oxo_cyclic)
    all_three_oxo = list(matches_standard) + list(matches_cyclic)
    
    if not all_three_oxo:
        return False, "3-oxo fatty acyl fragment not found"
    
    # Verify that at least one 3-oxo fragment is directly attached to the CoA moiety.
    # We expect that the sulfur atom in the 3-oxo fragment (last atom in our SMARTS pattern)
    # is connected to an atom from the CoA fragment.
    attached = False
    for match in all_three_oxo:
        # For standard pattern, the S is the fourth atom; for the cyclic one, use the matched S (position may be 0)
        # We check all atoms in the match: if any neighbor is in coa_atoms, assume proper connection.
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in coa_atoms:
                    attached = True
                    break
            if attached:
                break
        if attached:
            break

    if not attached:
        return False, "3-oxo acyl fragment found but not attached to the CoA moiety"

    return True, "Contains CoA moiety and 3-oxo fatty acyl fragment attached via thioester bond"

# Example usage:
if __name__ == "__main__":
    # Test one of the provided SMILES, e.g., (7Z)-3-oxohexadecenoyl-CoA
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)CC(CCC/C=C\\CCCCCCCC)=O)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_3_oxo_fatty_acyl_CoA(test_smiles)
    print(result, reason)