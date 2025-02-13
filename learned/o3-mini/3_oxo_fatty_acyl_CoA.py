"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: 3-oxo-fatty acyl-CoA

An oxo fatty acyl-CoA results from the condensation of the thiol group of coenzyme A
with the carboxyl group of any 3-oxo fatty acid. A key motif is that the fatty acyl
chain contains two carbonyl groups (one at the acid end and one at the 3-oxo substituent)
separated by one intervening carbon, and it is attached via a thioester bond to the CoA moiety.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    
    The method performs the following steps:
      1. Parse the SMILES and check that the overall formal charge is zero 
         (by summing the formal charges of all atoms).
      2. Search for a Coenzyme A fragment using a characteristic SMARTS pattern.
      3. Search for a 3-oxo fatty acyl fragment. This is defined by an open-chain motif
         –C(=O)[CX4]C(=O)[S]– (or a slightly different cyclic variant). The pattern should
         be directly attached (via a neighboring atom) to the CoA moiety.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        (bool, str):
            bool: True if the molecule is classified as a 3-oxo-fatty acyl-CoA; otherwise False.
            str: A brief explanation of the classification result.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Instead of mol.GetFormalCharge(), sum the formal charges of all atoms.
    overall_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if overall_charge != 0:
        return False, f"Molecule carries a non-zero overall formal charge ({overall_charge})"
    
    # Define a SMARTS pattern to detect the Coenzyme A moiety.
    # This pattern is heuristic and common to many CoA derivatives.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "Coenzyme A moiety not found"
    
    # Create a set of atom indices that are part of any CoA match.
    coa_atoms = set()
    for match in coa_matches:
        coa_atoms.update(match)
    
    # Define the standard open-chain 3-oxo fatty acyl fragment pattern:
    # It covers –C(=O)[CX4]C(=O)[S]– (the [CX4] atom represents an sp3 carbon, including potential branching).
    three_oxo_pattern = Chem.MolFromSmarts("C(=O)[CX4]C(=O)[S]")
    # Also, allow for a cyclic variant (e.g. a cyclohexanone derivative attached via thioester).
    three_oxo_cyclic = Chem.MolFromSmarts("[S]C(=O)C1=CCCCC1=O")
    
    # Retrieve any substructure matches for the above patterns.
    matches_standard = mol.GetSubstructMatches(three_oxo_pattern)
    matches_cyclic = mol.GetSubstructMatches(three_oxo_cyclic)
    all_three_oxo_matches = list(matches_standard) + list(matches_cyclic)
    
    if not all_three_oxo_matches:
        return False, "3-oxo fatty acyl fragment not found"
    
    # Verify that at least one 3-oxo fragment is directly attached to the CoA moiety.
    attached = False
    for match in all_three_oxo_matches:
        # Traverse each atom in the match.
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check if any neighbor of this atom is in the CoA fragment.
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
    # Test with one of the provided SMILES, e.g., (7Z)-3-oxohexadecenoyl-CoA.
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)CC(CCC/C=C\\CCCCCCCC)=O)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_3_oxo_fatty_acyl_CoA(test_smiles)
    print(result, reason)