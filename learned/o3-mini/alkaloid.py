"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: Alkaloid
Definition: Any of the naturally occurring, basic nitrogen compounds (mostly heterocyclic) found
mostly in the plant kingdom (but also in bacteria, fungi, and animals). By extension, certain neutral
compounds biogenetically related to basic alkaloids are also classed as alkaloids.
Note: Compounds in which the nitrogen is exocyclic (for example dopamine, mescaline, serotonin, etc.)
are usually classed as amines – however, some molecules (e.g. selegiline, ephedrine) are accepted alkaloids.
This heuristic now requires that at least one nitrogen is either (a) part of a ring or (b) 
exocyclic but directly attached to an aromatic atom.
"""
from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a given molecule (SMILES) is likely an alkaloid using an improved heuristic.
    The heuristic criteria are:
      1) SMILES must be valid (i.e. parseable).
      2) The molecule must have at least one nitrogen atom.
      3) At least one nitrogen atom must either be part of a ring (typically heterocyclic) OR, if exocyclic,
         it should be directly bonded to an aromatic atom.
         
    Args:
      smiles (str): SMILES string of the molecule
      
    Returns:
      bool: True if the molecule is alkaloid-like, False otherwise.
      str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    total_nitrogen = 0            # count all nitrogen atoms
    ring_nitrogen = 0             # nitrogen atoms that are in a ring
    exocyclic_aromatic_nitrogen = 0  # nitrogen not in a ring but bonded to at least one aromatic atom

    # Iterate over all atoms: count nitrogens and check ring membership or aromatic neighbor
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # nitrogen atomic number
            total_nitrogen += 1
            if atom.IsInRing():
                ring_nitrogen += 1
            else:
                # For exocyclic nitrogen, see if any neighbor is an aromatic atom.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIsAromatic():
                        exocyclic_aromatic_nitrogen += 1
                        break  # count each nitrogen only once

    if total_nitrogen == 0:
        return False, "No nitrogen atoms present, so unlikely to be an alkaloid"
      
    # Count as "alkaloid-like" if at least one nitrogen meets our criteria:
    qualifying_nitrogens = ring_nitrogen + exocyclic_aromatic_nitrogen
    
    if qualifying_nitrogens > 0:
        return True, f"Found {ring_nitrogen} ring nitrogen(s) and {exocyclic_aromatic_nitrogen} exocyclic nitrogen(s) attached to an aromatic system out of {total_nitrogen} nitrogen(s): consistent with an alkaloid classification"
    else:
        return False, "Nitrogen atoms are present but none appear in a ring or attached to an aromatic system; suggests a simple amine rather than an alkaloid"

# Example usage:
if __name__ == "__main__":
    # Try one of the examples from the provided list (e.g. selegiline)
    test_smiles = "[H][N+](C)(CC#C)C(C)Cc1ccccc1"  # selegiline(1+)
    result, reason = is_alkaloid(test_smiles)
    print("Alkaloid:", result)
    print("Reason:", reason)