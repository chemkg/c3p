"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide (gamma‐lactone with a 2‐furanone skeleton and substituted derivatives)

A butenolide (or 2‐furanone) is defined here as a five‐membered lactone ring that contains exactly one ring oxygen 
and a lactone (C=O) function built into the ring. 
This implementation searches for a substructure matching the canonical 2‐furanone core.
    
Requires: rdkit
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide (gamma-lactone with a 2-furanone skeleton)
    based on its SMILES string.

    The approach is simplified by trying to match the well‐known 2‐furanone substructure 
    using SMARTS. The SMARTS pattern "O=C1OC=CC1" represents a five‐membered ring where one 
    atom is an oxygen, one carbon bears the carbonyl, and two of the ring bonds are double. 
    (Substituents on the ring do not affect the match.)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule appears to be a butenolide, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the canonical 2-furanone (butenolide) skeleton.
    # This pattern represents a five-membered ring in which:
    #   – a carbon is double-bonded to oxygen (lactone carbonyl),
    #   – an oxygen is present in the ring,
    #   – and there is at least one C=C double bond within the ring.
    #
    # Note: The substructure match is tolerant of additional substituents attached to the ring atoms.
    butenolide_smarts = "O=C1OC=CC1"
    query = Chem.MolFromSmarts(butenolide_smarts)
    if query is None:
        return False, "Error in SMARTS pattern for butenolide"

    # Use RDKit substructure matching.
    if mol.HasSubstructMatch(query):
        return True, "Contains a 2-furanone substructure typical of butenolide derivatives"
    else:
        return False, "Does not contain the required 2-furanone skeleton for butenolides"

# Example usage (for testing purposes, remove or comment out in production)
if __name__ == '__main__':
    test_examples = {
        "protoanemonin": "C=C1OC(=O)C=C1",
        "Aspersclerotiorone B": "O=C1O/C(=C(/C[C@]2(OC=C(C2=O)C)C)\\C)/C(=C1)OC",
        "deslanoside": "[H][C@@]1(C[C@H](O)[C@]([H])(O[C@H]2C[C@H](O)[C@]([H])(O[C@H]3C[C@H](O)[C@]([H])(O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)[C@@H](C)O3)[C@@H](C)O2)[C@@H](C)O1)O[C@H]1CC[C@@]2(C)[C@]([H])(CC[C@]3([H])[C@]2([H])C[C@@H](O)[C@]2(C)[C@H](CC[C@]32O)C2=CC(=O)OC2",
        "5-Pentyl-3h-furan-2-one": "O1C(CCCCC)=CCC1=O",
        "Non-butenolide ester": "CC(=O)OC1=CC=CC=C1",
    }
    
    for name, smi in test_examples.items():
        result, reason = is_butenolide(smi)
        print(f"{name}: {result} – {reason}")