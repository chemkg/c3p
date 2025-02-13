"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: 3-oxo-fatty acyl-CoA

An oxo fatty acyl-CoA results from the condensation of the thiol group of
coenzyme A with the carboxyl group of a 3-oxo fatty acid. The hallmark motif
is that the fatty acyl chain contains two carbonyl groups separated by one 
intervening aliphatic carbon, and the acyl chain is attached via a thioester bond 
to the CoA moiety.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    
    The steps are:
      1. Parse the SMILES and check overall formal charge is zero.
      2. Look for a Coenzyme A moiety using a characteristic SMARTS pattern.
      3. Look for a 3-oxo fatty acyl fragment. Here we demand that a 
         linear (or near linear) acyl fragment matches a pattern:
           [#6;!A](=O) [C;!A] [C;!A](=O)[S]
         That is, a carbonyl carbon (non‐aromatic), a bridging aliphatic carbon,
         and a second carbonyl directly attached to a thioester S atom.
         For at least one such match we further require that one of its atoms is
         directly bonded to an atom from the CoA fragment.
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str):
          bool: True if the molecule is classified as a 3-oxo-fatty acyl-CoA,
                False otherwise.
          str: A brief explanation.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall charge (sum of formal charges over all atoms)
    overall_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if overall_charge != 0:
        return False, f"Molecule carries a non-zero overall formal charge ({overall_charge})"
    
    # Use a Coenzyme A fragment SMARTS.
    # This pattern is heuristic and common to many CoA derivatives.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "Coenzyme A moiety not found"
    
    # Gather set of atom indices that are part of any CoA match.
    coa_atoms = set()
    for match in coa_matches:
        coa_atoms.update(match)
    
    # Define a stricter SMARTS pattern for the 3-oxo fatty acyl fragment.
    # This pattern schematically looks for:
    #   aliphatic carbonyl carbon (non-aromatic) - bridging aliphatic carbon - 
    #   aliphatic carbonyl carbon (non-aromatic) attached to a sulfur.
    # The SMARTS below matches four atoms in order:
    #   [#6;!a](=O) - [C;!a] - [#6;!a](=O) - [S]
    three_oxo_pattern = Chem.MolFromSmarts("[#6;!a](=O)[C;!a][#6;!a](=O)[S]")
    
    # We also allow a cyclic variant in which the 3-oxo motif is embedded
    # in a ring. (This is a fallback if the standard pattern is not found.)
    three_oxo_cyclic = Chem.MolFromSmarts("[S][C;!a](=O)C1CCCC1=O")
    
    # Get all matches for both patterns.
    matches_standard = mol.GetSubstructMatches(three_oxo_pattern)
    matches_cyclic = mol.GetSubstructMatches(three_oxo_cyclic)
    all_matches = list(matches_standard) + list(matches_cyclic)
    
    if not all_matches:
        return False, "3-oxo fatty acyl fragment not found"
    
    # Now, for each match, we further check that:
    #   (a) None of the atoms in the match (the acyl fragment) is aromatic.
    #   (b) At least one atom from the match is directly bonded to a CoA atom.
    valid_match_found = False
    for match in all_matches:
        # Skip if any of the atoms in the match is aromatic.
        skip = False
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetIsAromatic():
                skip = True
                break
        if skip:
            continue
        
        # Check if any atom in this acyl fragment is directly attached
        # to an atom from the CoA fragment.
        attached = False
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in coa_atoms:
                    attached = True
                    break
            if attached:
                break
        
        if attached:
            valid_match_found = True
            break
    
    if not valid_match_found:
        return False, "3-oxo acyl fragment found but not directly attached to the CoA moiety or not a valid aliphatic fragment."
    
    return True, "Contains CoA moiety and valid 3-oxo fatty acyl fragment attached via thioester bond"

# Example usage:
if __name__ == "__main__":
    # Test with one SMILES (one of the true positives in the outcomes)
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)CC(CCC/C=C\\CCCCCCCC)=O)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_3_oxo_fatty_acyl_CoA(test_smiles)
    print(result, reason)