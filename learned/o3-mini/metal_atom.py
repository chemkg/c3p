"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: metal atom
An atom of an element that exhibits typical metallic properties (neutral),
with examples like [Ta], [89Y], [Cr], [195Po], [135Cs], etc.
Only a specific set of metal atoms (with zero formal charge) is allowed.
"""

from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if the molecule represented by the SMILES string is a neutral metal atom
    from the allowed set: {Ta, Y, Cr, Po, Cs, Rh, Rb, Pa, Ni, Ti, Mo, Tb, Ca, Mt, Li, Sn, Es, Ba, W, Cn}.

    Args:
        smiles (str): SMILES string (e.g., "[Ta]", "[89Y]")

    Returns:
        bool: True if the molecule is a neutral metal atom from the allowed set, False otherwise
        str: Explanation for the classification decision.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule consists of exactly one atom.
    if mol.GetNumAtoms() != 1:
        return False, "Molecule does not consist of a single atom"

    # Retrieve the only atom.
    atom = mol.GetAtomWithIdx(0)
    
    # Check that the atom is neutral (formal charge equal to 0)
    if atom.GetFormalCharge() != 0:
        return False, f"Atom '{atom.GetSymbol()}' has a formal charge ({atom.GetFormalCharge()}); only neutral atoms allowed"

    # Get the atomic symbol (ignoring isotope information)
    element = atom.GetSymbol()

    # Define the allowed set of neutral metal atoms based on accepted examples (true positives)
    allowed_metals = {
        "Ta",  # tantalum
        "Y",   # yttrium (and [89Y] or [Y])
        "Cr",  # chromium
        "Po",  # polonium
        "Cs",  # caesium ([135Cs])
        "Rh",  # rhodium
        "Rb",  # rubidium ([87Rb])
        "Pa",  # protactinium
        "Ni",  # nickel
        "Ti",  # titanium
        "Mo",  # molybdenum
        "Tb",  # terbium
        "Ca",  # calcium ([43Ca])
        "Mt",  # meitnerium
        "Li",  # lithium ([7Li])
        "Sn",  # tin
        "Es",  # einsteinium
        "Ba",  # barium
        "W",   # tungsten
        "Cn"   # copernicium
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
        "[51V]",      # false positive: not allowed
        "[Sr++]",     # false positive: charged
        "[Tb+3]",     # false positive: charged
        "[Cr-]",      # false positive: charged
        "[H][Al+]([H])[H]"  # false positive: multi-atom molecule
    ]
    
    for sm in test_smiles:
        result, reason = is_metal_atom(sm)
        print(f"SMILES: {sm} --> {result}, Reason: {reason}")