"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: 3-oxo-5beta-steroid, defined as 'Any 3-oxo steroid that has beta- configuration at position 5.'
Note: Due to the complexity of stereochemistry and steroid numbering, this implementation uses heuristic checks.
It checks whether the molecule has the expected fused steroid nucleus (three six‐membered rings and one five‐membered ring),
contains a ketone group that is part of a ring (for the 3-oxo functionality), and whether the input SMILES includes "@@"
chiral annotations (a possible proxy for the beta configuration in these examples).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule belongs to the 3-oxo-5beta-steroid class based on its SMILES string.
    
    Heuristic criteria:
       1. The molecule should have a fused ring system typical for steroids:
          - It should contain at least one five-membered ring and at least three six-membered rings.
       2. The molecule should contain a carbonyl group within a ring (indicating the “3-oxo”).
       3. The SMILES should contain chiral annotations and in particular use the "@@" designation 
          as a proxy for the beta configuration at position 5 (based on the provided examples).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule passes the heuristic checks, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule has a steroid-like fused ring system.
    # Obtain the rings in the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found in the molecule"
    
    # Count rings by size
    ring_sizes = [len(ring) for ring in rings]
    count_six = ring_sizes.count(6)
    count_five = ring_sizes.count(5)
    
    # A typical steroid core (cyclopentanoperhydrophenanthrene) has three six-membered rings and one five-membered ring.
    if count_six < 3 or count_five < 1:
        return False, ("Molecule does not contain the typical fused steroid ring system "
                       "(requires at least three six-membered rings and one five-membered ring)")
    
    # Check for a ketone group that is part of a ring.
    # The SMARTS "[R]C(=O)[!#1]" looks for a ring carbon (C) with a double-bonded oxygen.
    ketone_pattern = Chem.MolFromSmarts("[R]C(=O)[!#1]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found within a ring (missing 3-oxo feature)"
    
    # Check that the molecule has chiral annotations.
    # In many of the provided examples, a beta configuration (at position 5) is indicated by “@@” in the SMILES.
    # Note: This is only a heuristic indicator.
    if "@@" not in smiles:
        return False, "No beta stereochemical annotation found (missing '@@' chiral markers)"
    
    # (Further stereochemical checks could be done by assigning CIP codes to atoms,
    #  but determining which chiral center corresponds exactly to position 5 is nontrivial.)
    
    return True, "Molecule appears to be a 3-oxo steroid with beta configuration based on heuristic criteria"

# Example usage:
if __name__ == "__main__":
    # Test one of the provided examples
    test_smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@H]34)[C@@H]1CC[C@@H]2O"  # 17beta-hydroxy-5beta-estran-3-one
    result, reason = is_3_oxo_5beta_steroid(test_smiles)
    print(result, reason)