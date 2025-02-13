"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: Mineral nutrient
Definition: A mineral that is an inorganic nutrient which must be ingested and absorbed in
adequate amounts to satisfy a wide range of essential metabolic and/or structural functions 
in the human body.

The heuristic here is that mineral nutrients (at least among our examples) are salts or complexes 
that contain a metal (or metalloid) ion. We check if the molecule contains at least one atom 
whose elemental symbol is one of a set of common mineral nutrient metals/metalloids.
"""

from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines whether a molecule (given as SMILES) is a mineral nutrient.
    We use the heuristic that a mineral nutrient must contain at least one metal (or metalloid)
    ion from a predefined set.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        tuple: (bool, str) where bool indicates classification as a mineral nutrient and
               str gives the reason for the classification.
    """
    # Try to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a set of metal/metalloid ions that appear in our examples.
    # (This is a heuristic list based on the provided examples.)
    mineral_metals = {"Fe", "Zn", "Ca", "K", "Mg", "Na", "Ba", "Al", "Cs", "Pd", "Sb", "La"}

    # Scan molecule atoms for any metal from our list.
    found_metals = set()
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in mineral_metals:
            found_metals.add(symbol)
    
    # If at least one metal ion is found, we classify the compound as a mineral nutrient.
    if found_metals:
        return True, f"Contains mineral metal(s): {', '.join(sorted(found_metals))}"
    else:
        return False, "No recognized mineral metal found; not a mineral nutrient"
        
# For testing purposes, one might run e.g.:
# print(is_mineral_nutrient("[Fe+3].[O-]P([O-])(=O)[O-]"))
# print(is_mineral_nutrient("CC(=O)OC1=CC=CC=C1C(=O)O"))