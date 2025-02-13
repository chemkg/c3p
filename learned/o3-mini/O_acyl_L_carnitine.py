"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine 
Definition: An O-acylcarnitine in which the carnitine component has L-configuration.
We require that the molecule contains a carnitine substructure with the following features:
  - A chiral carbon (the carnitine center) which is attached to:
      • an oxygen that is esterified with a carbonyl group (the acyl substituent),
      • a carboxylate group (the propanoate part of carnitine),
      • and a trimethylammonium moiety.
  - When assigned CIP stereochemistry, the carnitine chiral center must have code “R” 
    (which corresponds to naturally occurring L-carnitine).
    
The improvements made compared to the earlier attempt are:
  • We use two SMARTS patterns (one for each chirality flag) but do not immediately return a negative
    result if one match has the wrong CIP code. We iterate over all matches and accept if any have the expected (R) configuration.
  • We define the SMARTS in a manner that does not require explicit hydrogen/deuterium annotations.
  
If no substructure is found, the function returns False with an explanation.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determine if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    The molecule must contain an O-acyl-carnitine skeleton and the carnitine chiral center must have
    the natural configuration for L-carnitine. In CIP naming this is "R".
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Assign stereochemistry so that CIP codes are set
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # The SMARTS below look for a carnitine-like fragment:
    #   • A chiral center (mapped as atom :1) that is attached to:
    #       - An oxygen that is further bonded to a carbonyl (O-acyl group).
    #       - A propanoate-like fragment: "CC(=O)[O-]" (the carboxylate group).
    #       - A trimethylammonium substituent: "C[N+](C)(C)C".
    #
    # We omit explicit isotopic specification so that deuterated forms still match.
    pattern_smarts1 = "[C@H:1]([O][C](=O)[*!#1])([C,N][N+](C)(C)C)(CC(=O)[O-])"
    pattern_smarts2 = "[C@@H:1]([O][C](=O)[*!#1])([C,N][N+](C)(C)C)(CC(=O)[O-])"
    # Note: The fragment [O][C](=O)[*!#1] means an oxygen bonded to a carbonyl carbon 
    # (the [*!#1] wildcard will match any non-hydrogen so that the acyl chain is free-form)
    
    patt1 = Chem.MolFromSmarts(pattern_smarts1)
    patt2 = Chem.MolFromSmarts(pattern_smarts2)
    
    if patt1 is None and patt2 is None:
        return False, "Could not compile SMARTS patterns"
    
    matches = []
    if patt1:
        matches.extend(mol.GetSubstructMatches(patt1, useChirality=True))
    if patt2:
        matches.extend(mol.GetSubstructMatches(patt2, useChirality=True))
    
    if not matches:
        return False, "No substructure matching an O-acyl-carnitine skeleton was found"
    
    # Instead of returning on the very first match, check all matches.
    valid_matches = 0
    wrong_config = 0
    for match in matches:
        # The first atom in the match (mapped from :1) is the carnitine center.
        try:
            chiral_atom = mol.GetAtomWithIdx(match[0])
        except IndexError:
            continue
        # Check if the CIP code property is available
        if not chiral_atom.HasProp("_CIPCode"):
            continue
        cip = chiral_atom.GetProp("_CIPCode")
        if cip == "R":
            valid_matches += 1
        else:
            wrong_config += 1
    
    if valid_matches > 0:
        return True, "Matches O-acyl-L-carnitine pattern with correct (R) stereochemistry at carnitine center"
    elif wrong_config > 0:
        return False, f"Found matches but the carnitine center(s) have wrong CIP configuration (none had 'R')"
    else:
        return False, "No match found with sufficient stereochemistry information"

# Example usage:
if __name__ == '__main__':
    # Test with one of the provided examples: O-acetyl-L-carnitine
    test_smiles = "CC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C"
    result, reason = is_O_acyl_L_carnitine(test_smiles)
    print("Result:", result)
    print("Reason:", reason)