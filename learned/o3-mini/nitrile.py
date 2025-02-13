"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: nitrile compounds (RC#N)
A nitrile is defined as a compound featuring a carbon–nitrogen triple bond 
where the carbon is substituted (i.e. not simply HC#N). The function disqualifies
nitrile groups that are coordinated to metals (e.g. in salt forms) but accepts any 
substituent (even heteroatoms) as long as it is not merely an implicit hydrogen.
"""

from rdkit import Chem

# Define a set of common metal symbols to filter out metal-bound nitriles.
metal_symbols = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac",
    "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
}

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile should contain at least one triple-bonded C#N group where the carbon is substituted.
    To be considered substituted, the nitrile carbon must have at least one additional explicit neighbor
    (i.e. a bond other than the one to the nitrile nitrogen). Also, if this additional neighbor is a metal,
    the group is not counted, as metals (or ionic compounds) are outside our target organic class.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a nitrile (contains at least one appropriate substituted nitrile group),
              False otherwise.
        str: Explanation of the reasoning and results.
    """
    # Parse the SMILES string, return error if invalid.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the basic nitrile SMARTS pattern (a carbon triple-bonded to nitrogen).
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if not nitrile_matches:
        return False, "No nitrile (C#N) functional group found"
    
    # Check each nitrile match.
    for match in nitrile_matches:
        # match[0] is the nitrile carbon; match[1] is the nitrile nitrogen.
        nitrile_c = mol.GetAtomWithIdx(match[0])
        nitrile_n = mol.GetAtomWithIdx(match[1])
        
        # Get the list of explicit neighbors of the nitrile carbon.
        # Note: RDKit does not count implicit hydrogens.
        neighbors = nitrile_c.GetNeighbors()
        # In an unsubstituted nitrile (HC#N), the nitrile carbon would have only one neighbor: the nitrile N.
        if len(neighbors) <= 1:
            # This nitrile is effectively HCN.
            continue
        
        # Identify the substituent neighbor(s) other than the nitrile nitrogen.
        # We permit any neighbor (even if heteroatom) as long as it is not a metal.
        valid_substituent_found = False
        for nb in neighbors:
            # Skip the nitrile nitrogen.
            if nb.GetIdx() == nitrile_n.GetIdx():
                continue
            # Check if the neighbor is a metal.
            if nb.GetSymbol() in metal_symbols:
                continue
            # If we find at least one non-metal neighbor (explicitly represented) then the nitrile carbon is substituted.
            valid_substituent_found = True
            break
        
        if valid_substituent_found:
            return True, "Molecule contains a substituted nitrile group (RC#N)"
    
    # If none of the nitrile matches met the requirement, then only unsubstituted nitrile or metal-bound nitrile was found.
    return False, "Only unsubstituted nitrile group(s) (HC#N) found; no C-substituted nitrile group present"

# Example usage (can be removed):
if __name__ == "__main__":
    test_smiles = "CC#N"  # Acetonitrile, a valid nitrile.
    result, reason = is_nitrile(test_smiles)
    print(test_smiles, result, reason)