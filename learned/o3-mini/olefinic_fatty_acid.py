"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Olefinic Fatty Acid
Definition: Any fatty acid (i.e. acyclic, long‐chain carboxylic acid) that contains at least one C=C double bond.
Additional criteria:
 - The carboxyl group should be terminal (its carbon has only one carbon neighbor).
 - The molecule should contain no rings.
 - The molecule should contain a minimal number (e.g. 8) of carbon atoms.
"""

from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    Olefinic fatty acids are defined as fatty acids (molecules containing a terminal carboxylic acid group,
    with an aliphatic long chain and no rings) that have at least one C=C double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is acyclic.
    # Fatty acids are long-chain aliphatic acids and should not contain rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a simple acyclic fatty acid"
    
    # Check for a carboxylic acid group.
    # This SMARTS pattern looks for a carbon that is double-bonded to an oxygen and
    # single-bonded to an -OH or its deprotonated equivalent.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # Verify that the carboxylic acid group is terminal.
    # In a fatty acid the acid carbon should have only one carbon neighbor.
    terminal_acid = False
    for match in acid_matches:
        acid_c_idx = match[0]  # first atom in the pattern is the acid carbon.
        # Retrieve neighboring atoms that are carbon.
        neighbors = [nbr for nbr in mol.GetAtomWithIdx(acid_c_idx).GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(neighbors) == 1:
            terminal_acid = True
            break
    if not terminal_acid:
        return False, "Carboxylic acid group is not terminal; not a typical fatty acid"
    
    # Check for a long aliphatic chain by counting carbon atoms.
    # We require at least 8 carbons.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 8:
        return False, "Too few carbon atoms to be a fatty acid"
    
    # Check for the presence of at least one carbon-carbon double bond.
    # The SMARTS "C=C" ensures that a C=C double bond is present
    # and it won't match the C=O since oxygen is not carbon.
    olefin_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(olefin_pattern):
        return False, "No C=C double bonds found; not an olefinic fatty acid"
    
    return True, "Contains a terminal carboxylic acid group, is acyclic with a long aliphatic chain, and has at least one C=C double bond, classifying it as an olefinic fatty acid"