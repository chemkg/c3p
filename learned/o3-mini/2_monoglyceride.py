"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride, defined as 'A monoglyceride in which the acyl substituent is located at position 2.'
The key fragment is a glycerol backbone where the central carbon (CH) is bonded to:
    • two CH2OH groups (the terminal hydroxyl groups) and
    • one O–C(=O)[*] group (the acyl substituent via an ester bond).
This program uses a SMARTS pattern to search for that characteristic fragment.
"""

from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has a glycerol backbone where the central (2-) carbon bears an acyl ester group.
    The key structural fragment is:
        CH2OH – CH(OC(=O)R) – CH2OH
    We define a SMARTS pattern that exactly captures a carbon attached to two CH2OH groups
    and one O–C(=O)[*] substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is identified as a 2-monoglyceride, False otherwise.
        str: Explanation for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the 2-monoglyceride substructure.
    # Explanation:
    #   [C]                -- any carbon (ignores stereochemical markers)
    #   ([CH2]O)          -- a CH2 group attached to an -OH (terminal CH2OH group)
    #   ([CH2]O)          -- a second terminal CH2OH group
    #   O[C](=O)[*]       -- an oxygen connected to a carbonyl carbon (i.e. ester bond) and any acyl chain
    # The entire pattern is enclosed in $([...]) so that it is treated as one unit/substructure.
    # This pattern should match exactly one occurrence if the molecule is a monoglyceride esterified at the 2-position.
    pattern_smarts = "[$([C]([CH2]O)([CH2]O)O[C](=O)[*])]"
    pattern = Chem.MolFromSmarts(pattern_smarts)
    if pattern is None:
        return False, "SMARTS pattern could not be compiled"
    
    # Find all matching substructure fragments.
    matches = mol.GetSubstructMatches(pattern)
    
    # We expect exactly one match for a proper 2-monoglyceride.
    if len(matches) == 0:
        return False, "2-monoglyceride substructure not found"
    elif len(matches) > 1:
        return False, "Multiple 2-monoglyceride-like substructures found; molecule may be di-/tri-acylated or not a simple 2-monoglyceride"
    
    # If exactly one match is found, we classify the molecule as a 2-monoglyceride.
    return True, "Molecule contains a glycerol backbone with a single acyl ester at the 2-position"

# For testing purposes (you can remove or comment these lines in production):
if __name__ == "__main__":
    test_smiles = [
        "O=C(CCC/C=C\\C\\C=C/C=C/C(C/C=C\\CCCCC)OO)OC(CO)CO",  # 12-HPETE 2-glyceryl ester
        "CCCCCCCCCCCCCCCC(=O)OC(CO)CO",  # 2-palmitoylglycerol
        "O(C(=O)CCCCCCCCCCCCCCCCCCCCC)C(CO)CO",  # MG(0:0/22:0/0:0)
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP([O-])(=O)OC[C@H]([NH3+])C([O-])=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC" # false positive example
    ]
    for s in test_smiles:
        flag, reason = is_2_monoglyceride(s)
        print(f"SMILES: {s}\nResult: {flag}\nReason: {reason}\n")