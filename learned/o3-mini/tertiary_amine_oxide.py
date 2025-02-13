"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: Tertiary Amine Oxide
A tertiary amine oxide is defined as an N-oxide where there are three organic groups bonded to the nitrogen atom.
This means that in the molecule there should be a nitrogen atom that is bonded to one oxygen (having formal charge -1)
and three other substituents. The nitrogen atom itself carries a formal positive charge.
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide contains an oxidized, positively charged nitrogen bonded to one oxygen (negatively charged)
    and three organic substituents (we require these substituents to be directly bonded to a carbon atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a tertiary amine oxide, otherwise False.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Loop over all atoms to look for an N-oxide that meets tertiary amine criteria.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        # Check that nitrogen carries a positive charge.
        if atom.GetFormalCharge() != 1:
            continue
        # A tertiary amine oxide nitrogen is bonded to 4 atoms (one O- and three organic groups).
        if atom.GetDegree() != 4:
            continue

        neighbors = list(atom.GetNeighbors())
        oxide_found = False
        organic_neighbors = []
        # Loop over neighbors and sort into oxygen (that should be O-) and others.
        for neigh in neighbors:
            if neigh.GetSymbol() == "O" and neigh.GetFormalCharge() == -1:
                oxide_found = True
            else:
                organic_neighbors.append(neigh)

        if not oxide_found:
            continue
        # There must be exactly three organic substituents (other than the oxygen).
        if len(organic_neighbors) != 3:
            continue

        # Check that each organic substituent is carbon-based.
        # For simplicity we require that the directly bonded neighbor is a carbon atom.
        for sub in organic_neighbors:
            if sub.GetSymbol() != "C":
                return False, "Found an N-oxide but one substituent is not carbon-based"
        # If we get here then we have located a tertiary amine oxide group.
        return True, "Found tertiary amine oxide with N(+)-O(-) and three carbon containing substituents."

    # No tertiary amine oxide group found.
    return False, "No tertiary amine oxide pattern found"

# Example usage (uncomment to test):
# test_smiles = "C[N+](C)([O-])C"  # trimethylamine N-oxide
# result, reason = is_tertiary_amine_oxide(test_smiles)
# print(result, reason)