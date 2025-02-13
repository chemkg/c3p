"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: metal atom
An atom of an element that exhibits typical metallic properties (neutral and unbonded).
Examples include individual atoms such as [Ta], [89Y], [Cr], [195Po], [135Cs], etc.
This classifier only accepts single, unbonded, neutral atoms from a standard set of metals.
"""

from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if the molecule represented by the SMILES string is a neutral, unbonded metal atom.
    
    Criteria:
      1. The molecule must consist of exactly one atom (and no bonds).
      2. The atom must have a formal charge of 0.
      3. Its elemental symbol (ignoring any isotopic information) is in the allowed set of metal atoms.
      
    The allowed metal atoms are defined using a complete set of metals from the periodic table that
    generally exhibit metallic properties.
      
    Args:
        smiles (str): SMILES string (e.g., "[Ta]", "[89Y]")
        
    Returns:
        bool: True if the input is a neutral, single-atom metal from the allowed set, else False.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule consists of exactly one atom and has no bonds.
    if mol.GetNumAtoms() != 1 or mol.GetNumBonds() != 0:
        return False, "Molecule is not a single, unbonded atom"
    
    # Retrieve the only atom
    atom = mol.GetAtomWithIdx(0)
    
    # Check that the atom is neutral (formal charge 0)
    if atom.GetFormalCharge() != 0:
        return False, f"Atom '{atom.GetSymbol()}' has a formal charge ({atom.GetFormalCharge()}); only neutral atoms allowed"
    
    # Get the element symbol (ignoring isotope information)
    element = atom.GetSymbol()
    
    # Define a more complete allowed set of metal atoms.
    allowed_metals = {
        # Alkali metals
        "Li", "Na", "K", "Rb", "Cs", "Fr",
        # Alkaline earth metals
        "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
        # Transition metals
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        # Post-transition metals
        "Al", "Ga", "In", "Sn", "Tl", "Pb", "Bi",
        # Lanthanides
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        # Actinides
        "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
        # Special cases (other metals mentioned in examples)
        "Po", "Mt", "Cn"
    }
    
    if element not in allowed_metals:
        return False, f"Atom '{element}' is not in the allowed set of metal atoms"
    
    return True, f"Atom '{element}' is classified as a metal"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = [
        "[Ta]",
        "[89Y]",
        "[Cr]",
        "[195Po]",
        "[135Cs]",
        "[Rh]",
        "[216Po]",
        "[87Rb]",
        "[Pa]",
        "[Ni]",
        "[Ti]",
        "[95Mo]",
        "[Tb]",
        "[217Po]",
        "[43Ca]",
        "[Mt]",
        "[7Li]",
        "[116Sn]",
        "[194Po]",
        "[Es]",
        "[Y]",
        "[Ba]",
        "[201Po]",
        "[W]",
        "[Cn]",
        "[Yb]",   # now allowed
        "[Lu]",   # now allowed
        "[39K]",  # now allowed (potassium)
        "[Ra]",   # now allowed
        "[139La]",# now allowed (lanthanum)
        "[55Mn]", # now allowed (manganese)
        "[59Co]", # now allowed (cobalt)
        "[Bk]",   # now allowed (berkelium)
        "[Cm]",   # now allowed (curium)
        "[Al]",   # now allowed (aluminium)
        "[Li][H+]",         # false positive: multi-atom molecule
        "[H][Sn]([H])[H]",   # false positive: multi-atom molecule
        "[H][Po][H]",        # false positive: multi-atom molecule
        "[Li][H]",           # false positive: multi-atom molecule
        "[Li][H-]"           # false positive: multi-atom molecule
    ]
    
    for sm in test_smiles:
        result, reason = is_metal_atom(sm)
        print(f"SMILES: {sm} --> {result}, Reason: {reason}")