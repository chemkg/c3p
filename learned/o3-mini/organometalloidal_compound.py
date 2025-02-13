"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: Organometalloidal Compound
Definition: A compound having bonds between one or more metalloid (here, specifically arsenic) atoms and one or more carbon atoms of an organyl group.
We improve upon a simple check by only considering arsenic (As, atomic number 33) and by ensuring that if an As atom is attached to more than one carbon, at least one of the carbon substituents is larger than a simple methyl group.
"""

from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    For our implementation, we focus on arsenic-containing compounds. The molecule is classified as an organometalloidal compound if it contains at least one bond
    between an arsenic atom and a carbon atom, and if the carbon (or at least one of them, when there are multiple) can be considered an organyl group.
    
    An organyl group is here defined as either a solitary carbon substituent (e.g. a methyl group) OR, if an arsenic is attached to more than one carbon, at least one of the carbons 
    should not be "methyl" (i.e. it should have more than one heavy-atom neighbor or be aromatic).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule qualifies as an organometalloidal compound, False otherwise.
        str: Reason for the classification.
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # For our improved classification, only consider arsenic (atomic number 33)
    AS_ATOMIC_NUM = 33

    # Helper function to decide if a carbon atom is "methyl"
    def is_methyl(carbon):
        # A carbon will be considered a methyl if it is aliphatic (non-aromatic)
        # and its heavy-atom degree (neighbors excluding hydrogens) is exactly 1.
        # Note: In methylarsonic acid the only C is CH3 and that should count.
        if not carbon.GetSymbol() == "C":
            return False
        # Count heavy neighbors
        heavy_neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # Aromatic carbons we'll treat as non-methyl (being part of a larger aromatic system).
        if carbon.GetIsAromatic():
            return False
        return len(heavy_neighbors) == 1

    # Flag to indicate we found a valid As-C bond with proper organyl group
    found_valid_bond = False

    # Iterate over all atoms; focus on arsenic atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == AS_ATOMIC_NUM:
            # Get all carbon neighbors of this arsenic
            carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if not carbon_neighbors:
                continue  # no C connected to this As

            # If there is exactly one attached carbon, we accept it (even if it is a methyl)
            if len(carbon_neighbors) == 1:
                neighbor = carbon_neighbors[0]
                found_valid_bond = True
                return True, f"Found a bond between As and C ({neighbor.GetSmarts() if hasattr(neighbor, 'GetSmarts') else neighbor.GetSymbol()})."
            else:
                # More than one C is attached.
                # We require that at least one of these carbons is NOT just a methyl.
                for neighbor in carbon_neighbors:
                    if not is_methyl(neighbor):
                        found_valid_bond = True
                        return True, f"Found a bond between As and C where the C is part of a larger organyl group: {neighbor.GetSymbol()}."
                # If all attached carbons are simply methyl groups (and more than one is present) then we consider it not an organometalloidal compound.
                # This helps to weed out cases like di­methylarsinate.
                return False, "Arsenic is bonded to multiple only-methyl groups; does not qualify as an organyl bond."

    if not found_valid_bond:
        return False, "No bond found between arsenic and a carbon of an organyl group."
    
    # Fallback
    return False, "No valid organometalloidal (As–C) bond detected."

# Example usage:
if __name__ == "__main__":
    test_examples = [
        # True positives (should qualify):
        ("C[As](O)(O)=O", "methylarsonic acid"),
        ("[As](=O)(CCCCCCCCCCCCCCCCC)(C)C", "1-dimethylarsinoyl-heptadecane"),
        ("OC(=O)C[As](O)(O)=O", "arsenoacetic acid"),
        # False positive example: dimethylarsinate should not qualify.
        ("C[As](C)([O-])=O", "dimethylarsinate"),
        # False positive example: organosilicon compound should not qualify.
        ("CC[Si](C)(C)C", "ethyl(trimethyl)silane")
    ]
    for smi, name in test_examples:
        res, reason = is_organometalloidal_compound(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nRESULT: {res}\nREASON: {reason}\n{'-'*60}")