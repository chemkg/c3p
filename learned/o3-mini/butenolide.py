"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide (gamma-lactone with a 2-furanone skeleton and substituted derivatives)
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is defined as a gamma-lactone that consists of a 2-furanone skeleton
    (i.e. a five-membered cyclic ester ring with one ring oxygen, a carbonyl group, 
    and an internal double bond) and its substituted derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a butenolide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a basic 2-furanone skeleton:
    # This pattern looks for a five-membered ring "C1OC(=O)C=C1" that represents
    # a furanone core. It requires an oxygen in the ring, a carbonyl group, and a C=C bond.
    #
    # Note: Many butenolide structures are substituted. This pattern will match the 
    # canonical core motif. If your molecules employ variants that change the ring bonding,
    # you might wish to adapt this SMARTS.
    butenolide_pattern = Chem.MolFromSmarts("C1OC(=O)C=C1")
    
    # Search for the butenolide substructure in the molecule
    if mol.HasSubstructMatch(butenolide_pattern):
        return True, "Contains a 2-furanone skeleton typical of butenolide derivatives"
    else:
        return False, "Does not contain the required 2-furanone skeleton for butenolides"
        
# Example usage (you can remove or comment out these examples for production):
if __name__ == '__main__':
    examples = {
        "protoanemonin": "C=C1OC(=O)C=C1",
        "Appenolide A": "O=C1OCC(=C1C)/C=C/CCCCCC(=O)C",
        "Non-butenolide": "CC(=O)OC1=CC=CC=C1"  # example ester that lacks furanone ring
    }
    
    for name, smi in examples.items():
        result, reason = is_butenolide(smi)
        print(f"{name}: {result} – {reason}")