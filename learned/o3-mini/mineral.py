"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: minerals
A mineral is normally an inorganic crystalline (or formerly crystalline) substance,
often formed by geological processes. Heuristics used:
  - Many minerals are composed of ionic fragments (multiple disconnected components).
  - Many minerals include metal atoms.
  - Many show formal charges.
  - Many are relatively “inorganic” (e.g. lacking hydrogens or carbon–hydrogen bonds)
  
Note: Because crystallinity/geological history is not encoded in SMILES,
this classifier uses structural proxies and may not cover all edge cases.
"""

from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is likely a mineral based on its SMILES string using
    several heuristic checks (ionic nature, presence of metals, formal charges, etc.).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a mineral, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Try to parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get fragments: many ionic compounds are disconnected (e.g., salts)
    fragments = Chem.GetMolFrags(mol, asMols=True)
    multiple_fragments = len(fragments) > 1

    # Define a set of metal element symbols that are common in minerals.
    metal_symbols = {
        "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
        "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
        "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce",
        "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
        "Po", "Fr", "Ra", "Ac", "Th", "U", "Np", "Pu"
    }
    
    # Initialize flags for each heuristic condition
    metal_found = False
    has_formal_charge = False
    count_C = 0
    count_H = 0
    
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol == "C":
            count_C += 1
        if symbol == "H":
            count_H += 1
        # Check for metal atoms
        if symbol in metal_symbols:
            metal_found = True
        # Check if any atom has a non-zero formal charge
        if atom.GetFormalCharge() != 0:
            has_formal_charge = True

    # We also check for an "inorganic" signature:
    # Many minerals have either no carbon at all or very few hydrogens,
    # especially in salts such as carbonates or sulfates.
    no_carbon = (count_C == 0)
    no_hydrogen = (count_H == 0)
    
    # Determine if at least one of the conditions typical for a mineral is met.
    if metal_found:
        return True, "Contains metal element(s) typical of minerals."
    if multiple_fragments:
        return True, "Contains multiple disconnected fragments, consistent with an ionic compound."
    if has_formal_charge:
        return True, "Molecule contains formal charges, suggesting an ionic nature."
    if no_carbon:
        return True, "Molecule contains no carbon, a hallmark of many inorganic minerals."
    if no_hydrogen:
        return True, "Molecule contains no hydrogen, which is common for inorganic salts."
    
    # If none of our heuristics apply, then we are not classifying this as a mineral.
    return False, "Does not exhibit typical inorganic mineral characteristics (metal ions, ionic fragments, or absence of C/H)."
    
# Example usage:
if __name__ == '__main__':
    test_smiles = "[Fe++].[S-][S-]"  # pyrite
    is_min, reason = is_mineral(test_smiles)
    print("is_mineral:", is_min)
    print("Reason:", reason)