"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: Organic sulfide (thioether) compounds
Definition: Compounds having the structure R–S–R (where R ≠ H)
"""

from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.
    An organic sulfide is defined as a compound containing an R–S–R motif with both R groups
    being non-hydrogen substituents (i.e. no S–H bond is allowed).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as an organic sulfide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so we can check for S–H bonds.
    mol = Chem.AddHs(mol)
    
    # Iterate over all atoms to find sulfur atoms (atomic number 16).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:
            # For organic sulfide (thioether), the S should be bonded only to non-hydrogen atoms.
            # Count the number of heavy (non-hydrogen) neighbors.
            heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
            # We require at least 2 non-hydrogen bonds for an R-S-R motif.
            if len(heavy_neighbors) == 2:
                return True, "Contains an R-S-R motif with no S-H bonds (organic sulfide)"
    
    # If no S atom with exactly two non-hydrogen neighbors was found, then the molecule
    # does not contain the organic sulfide (thioether) signature.
    return False, "No R-S-R motif (thioether) found"

# Example test cases (uncomment to run tests)
# test_smiles = [
#     "S(CCCCC)C",  # 1-(Methylthio)pentane, should be True
#     "CSC",        # Dimethyl sulfide, should be True
#     "SC",         # Methanethiol, should be False since one substituent is H
#     "c1ccsc1"     # Thiophene, should be True (S in an aromatic ring bonded to 2 carbons)
# ]
#
# for s in test_smiles:
#     result, reason = is_organic_sulfide(s)
#     print(f"SMILES: {s} -> {result}: {reason}")