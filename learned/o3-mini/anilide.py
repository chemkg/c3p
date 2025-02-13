"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: Anilide 
Definition: Anilide is any aromatic amide obtained by acylation of aniline
Meaning that the amide (C(=O)-N) group must have the N attached to an aromatic (benzene-type) carbon.
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    Anilides are aromatic amides formed by acylation of aniline (i.e. the amide N must be directly bound to an aromatic ring).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an anilide, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use a SMARTS pattern to identify the amide bond.
    # The pattern [#6](=O)[#7] corresponds to a carbonyl carbon (atomic number 6) double-bonded to oxygen (implicitly)
    # and single-bonded to a nitrogen (atomic number 7). This will match any amide bond.
    amide_pattern = Chem.MolFromSmarts("[#6](=O)[#7]")
    matches = mol.GetSubstructMatches(amide_pattern)
    if not matches:
        return False, "No amide bond (C(=O)N) found"

    # For each amide match, check if the nitrogen (third atom in the SMARTS) is bound to an aromatic carbon.
    # This is to ensure the amide came from the acylation of aniline.
    for match in matches:
        # match order: (carbonyl carbon, oxygen, nitrogen)
        carbonyl_c = match[0]
        nitrogen = match[2]
        n_atom = mol.GetAtomWithIdx(nitrogen)
        # Check each neighbor of the amide nitrogen except the carbonyl carbon
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetIdx() == carbonyl_c:
                continue  # skip the carbonyl carbon, which is not the aromatic group from aniline
            # Check if the neighbor is aromatic (and a carbon)
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                return True, "Found aromatic amide bond indicative of an anilide structure"
    
    return False, "No amide bond with nitrogen attached to an aromatic ring was found"