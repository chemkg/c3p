"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester has the general structure -[N]-C(=O)-O-, where nitrogen may have substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for carbamate ester group, excluding ring nitrogens explicitly
    carbamate_pattern = Chem.MolFromSmarts("[N;!R]-C(=O)-[OX2]")
    # SMARTS for carbamate where the N is part of a ring structure
    carbamate_ring_pattern = Chem.MolFromSmarts("[N;R]-C(=O)-[OX2]")
    # SMARTS pattern for amides, to be excluded
    amide_pattern = Chem.MolFromSmarts("[N]-C(=O)")

    # Check for amides
    if mol.HasSubstructMatch(amide_pattern):
        for match in mol.GetSubstructMatches(amide_pattern):
            amide_atom = mol.GetAtomWithIdx(match[0])
            for nbr in amide_atom.GetNeighbors():
                if nbr.GetSymbol() == 'C':
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetSymbol() == 'O' and nbr2.GetIdx() in match:
                            
                            return False, "Contains an amide, not a carbamate"
    
    #Check for the carbamate pattern
    if mol.HasSubstructMatch(carbamate_pattern) or mol.HasSubstructMatch(carbamate_ring_pattern):
         return True, "Contains a carbamate ester group"
    else:
        return False, "Does not contain a carbamate ester group"