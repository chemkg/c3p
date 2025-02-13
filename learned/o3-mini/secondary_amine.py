"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: Secondary Amine
Definition: A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
This means that the nitrogen (N) should be bonded to exactly two carbon atoms (hydrocarbyl groups)
and have one hydrogen (either implicit or explicit).
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine features a nitrogen atom that has exactly one hydrogen and exactly two 
    carbon (hydrocarbyl) substituents.

    Args:
        smiles (str): SMILES string representing the molecule

    Returns:
        bool: True if the molecule contains at least one secondary amine functional group, False otherwise.
        str: A brief explanation of the classification decision.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms in the molecule looking for nitrogen atoms that might be secondary amines.
    found = False
    details = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'N':
            continue
        
        # Get the total number of hydrogens (both implicit and explicit) on the nitrogen atom.
        num_hydrogens = atom.GetTotalNumHs()
        
        # Count how many of the nitrogen's neighbors are carbon atoms.
        carbon_neighbors = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
        
        # A secondary amine (derived from ammonia by replacing two hydrogens) should have exactly:
        # - One hydrogen (num_hydrogens == 1) 
        # - Two carbon neighbors (carbon_neighbors == 2)
        if num_hydrogens == 1 and carbon_neighbors == 2:
            found = True
            details.append(f"Atom index {atom.GetIdx()} (N) has 1 hydrogen and 2 carbon substituents")
    
    if found:
        return True, "Secondary amine group found: " + "; ".join(details)
    else:
        return False, "No nitrogen with exactly one hydrogen and two carbon substituents (secondary amine) was found"
        
# Example usage:
if __name__ == "__main__":
    # Test examples: N-methylaniline and dimethylamine (which are secondary amines)
    test_smiles = [
        "CNc1ccccc1",   # N-methylaniline
        "[H]N(C)C",     # dimethylamine
        "CCNCC",        # diethylamine (tertiary amine) - should be false
        "c1cc[nH]c1"    # 1H-pyrrole (not an aliphatic secondary amine)
    ]
    
    for s in test_smiles:
        result, reason = is_secondary_amine(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")