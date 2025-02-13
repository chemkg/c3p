"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: Aldehyde (RC(=O)H)
Definition: A compound in which a carbonyl group (C=O) is bonded to exactly one hydrogen atom and one R group (which is not oxygen).
Note: Formaldehyde (H2C=O) and groups that are part of carboxylic acids (where the noncarbonyl substituent is oxygen) are excluded.
"""

from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is defined as a molecule containing a formyl group: a carbon atom that is 
    double-bonded to oxygen and that has exactly one hydrogen and one heavy (nonoxygen) substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an aldehyde, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Add explicit hydrogens so that the hydrogen count for the aldehyde carbon is accurate.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for the aldehyde group.
    # [CH1] ensures the carbon has exactly 1 hydrogen.
    # (=O) requires a double-bonded oxygen.
    # [!O] requires the other substituent to not be oxygen (thus avoiding acids, esters, etc.).
    aldehyde_smarts = "[CH1](=O)[!O]"
    pat = Chem.MolFromSmarts(aldehyde_smarts)
    
    # Check for a substructure match. This returns True if at least one aldehyde group is found.
    if mol.HasSubstructMatch(pat):
        return True, "Aldehyde group detected: a carbonyl carbon (C=O) bonded to exactly one hydrogen and one non-oxygen substituent."
    
    # If no match is found, then the molecule does not contain an aldehyde group as defined.
    return False, "No aldehyde functional group (RC(=O)H) detected."

# Example usage: testing some of the provided SMILES strings.
if __name__ == '__main__':
    smiles_examples = [
        "O=CC1=C(O)C=C(C#CC(=C)C)C(=C1)O",         # 2,5-dihydroxy-4-(3-methylbut-3-en-1-ynyl)benzaldehyde
        "C[C@H](C=O)C(O)=O",                         # (2R)-2-methyl-3-oxopropanoic acid (should be rejected: carboxylic acid group present)
        "CN(CCC=O)N=O",                              # 3-(N-Nitrosomethylamino)propionaldehyde
        "CC(C)[C@H](NC(=O)OCc1ccccc1)C(=O)N[C@@H](Cc1ccccc1)C=O",  # Z-Val-Phe-H
        "O=CCCCCCCCCC/C=C/C=C\\CC",                   # 11E,13Z-Hexadecadienal
        "[H]C(=O)CC#C",                              # but-3-ynal
        "[H]C(=O)C(C)CCCCCCCCCCCCC",                  # 2-methylpentadecanal
        "[H]C(=O)c1ccc2ccccc2c1",                     # 2-naphthaldehyde
        "[H]C(=O)C(C)C",                             # Simple aldehyde example: isobutyraldehyde (RC(=O)H)
        "OC[C@H](CC(C(=O)[H])=O)O"                    # (S)-4,5-dihydroxy-2-oxopentanal
    ]
    
    for s in smiles_examples:
        result, reason = is_aldehyde(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*60}")