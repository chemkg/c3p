"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies chemical entities of the class carbamate ester.
Definition: any ester of carbamic acid or its N‐substituted derivatives.
Characteristic substructure: R–O–C(=O)–NR′.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is defined as an ester of carbamic acid or its N-substituted derivatives,
    with the key substructure R-O-C(=O)-NR' (the O-C(=O)-N motif).
    
    This version uses a SMARTS pattern for [#8]-[#6](=O)-[#7] and then further checks that:
      - The oxygen (supposedly the ester oxygen) is esterified (i.e. has at least one heavy-atom neighbor besides the carbonyl carbon).
      - The carbonyl carbon is attached to exactly three atoms (the ester oxygen, the nitrogen, and the carbonyl oxygen via double bond).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a carbamate ester group, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern: oxygen bonded to a carbonyl carbon which in turn bonds to a nitrogen.
    # The pattern [#8]-[#6](=O)-[#7] matches any oxygen (atomic number 8), 
    # any carbon (atomic number 6) with a double-bonded oxygen, and any nitrogen (atomic number 7).
    pattern_smarts = "[#8]-[#6](=O)-[#7]"
    carbamate_pattern = Chem.MolFromSmarts(pattern_smarts)
    if carbamate_pattern is None:
        return False, "Error creating SMARTS pattern for carbamate ester"
    
    # Find substructure matches in the molecule.
    matches = mol.GetSubstructMatches(carbamate_pattern)
    if not matches:
        return False, "Carbamate ester group (O-C(=O)-N) not found in the molecule."
    
    # Loop over each match and perform additional checks for proper ester connectivity.
    for match in matches:
        # According to our SMARTS, match indices correspond to:
        #   match[0]: the oxygen atom (supposed ester oxygen)
        #   match[1]: the carbonyl carbon
        #   match[2]: the nitrogen atom
        o_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        n_atom = mol.GetAtomWithIdx(match[2])
        
        # Check 1: The oxygen should be esterified.
        # Typically an -OR group will have a degree of 2 or more (it is connected to one carbon and one additional heavy atom)
        # whereas a hydroxyl (-OH) typically has only one bond in the mol graph (the H is implicit).
        if o_atom.GetDegree() < 2:
            # This oxygen likely comes from a carboxylic acid rather than an ester.
            continue

        # Additionally, ensure that one of the oxygen neighbors (other than the acyl carbon) is a heavy atom and preferably a carbon.
        o_neighbors = [nbr for nbr in o_atom.GetNeighbors() if nbr.GetIdx() != c_atom.GetIdx()]
        if not any(nbr.GetAtomicNum() == 6 for nbr in o_neighbors):
            # If no neighbor is carbon, then the oxygen may not represent a typical ester linkage.
            continue
        
        # Check 2: The carbonyl carbon should have three connections:
        # It is double bonded to an oxygen (carbonyl oxygen) and single bonded to the ester oxygen and the nitrogen.
        # We do not count the double bond as two separate connections when using GetDegree (it counts as one bond per neighbor).
        if c_atom.GetDegree() != 3:
            continue
        
        # If both checks pass we assume a valid carbamate ester substructure is found.
        return True, "Molecule contains a carbamate ester group (O-C(=O)-N) with proper ester connectivity."
    
    # If none of the matches pass the additional criteria, then we consider the molecule not to contain a carbamate ester.
    return False, "Carbamate ester group (O-C(=O)-N) was found but none passed the connectivity checks for an ester."