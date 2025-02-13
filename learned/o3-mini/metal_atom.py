"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: metal atom
An atom of an element that exhibits typical metallic properties, being typically shiny,
with high electrical and thermal conductivity.
Examples: [Ta], [89Y], [Cr], [Po], [135Cs], etc.
"""

from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a molecule represented by a SMILES string is a metal atom.
    A metal atom is defined as a single atom (optionally with isotope notation) whose element is considered a metal.

    Args:
        smiles (str): SMILES string of the molecule (e.g., "[Ta]", "[89Y]")

    Returns:
        bool: True if the molecule is a metal atom, False otherwise.
        str: Reason for classification.
    """
    
    # Attempt to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule consists of exactly one atom
    if mol.GetNumAtoms() != 1:
        return False, "More than one atom present; not a single metal atom"

    # Retrieve the single atom from the molecule
    atom = mol.GetAtomWithIdx(0)
    element = atom.GetSymbol()  # Get the element symbol (isotope info is in atom.GetIsotope(), but symbol remains)
    
    # Define a set of element symbols that are considered metals.
    # This set includes alkali metals, alkaline earth metals, transition metals, post-transition metals,
    # lanthanides, and actinides (with polonium included).
    metal_elements = {
        # Alkali metals
        "Li", "Na", "K", "Rb", "Cs", "Fr",
        # Alkaline earth metals
        "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
        # Transition metals
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
        # Post-transition metals / metalloids that show metallic character
        "Al", "Ga", "In", "Sn", "Tl", "Pb", "Bi", "Po",
        # Lanthanides
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
        "Ho", "Er", "Tm", "Yb", "Lu",
        # Actinides
        "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
    }
    
    # Check if the atom's element is in the set of metals.
    if element not in metal_elements:
        return False, f"Atom '{element}' is not classified as a metal"
    
    return True, f"Atom '{element}' is classified as a metal"

# Example usage:
if __name__ == "__main__":
    test_smiles = ["[Ta]", "[89Y]", "[Cr]", "[Po]", "[135Cs]", "CCO"]  # last one should fail
    for sm in test_smiles:
        result, reason = is_metal_atom(sm)
        print(f"SMILES: {sm} --> {result}, Reason: {reason}")