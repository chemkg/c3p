"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: nitrile compounds (RC#N)
A nitrile is defined as a compound featuring a carbon–nitrogen triple bond 
where the carbon is substituted (i.e. not simply HC#N). This version converts 
all implicit hydrogens into explicit ones so that the degree of the nitrile carbon 
is counted correctly. A nitrile group is accepted only if on the nitrile carbon, 
besides the nitrile nitrogen, there is a neighbor that is neither hydrogen nor a metal.
"""

from rdkit import Chem

# Set of metal symbols that, if attached to a nitrile carbon, cause us to reject that group.
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
    Determines if a molecule is a substituted nitrile (RC#N) based on its SMILES string.
    
    The method adds explicit hydrogens so we can reliably count bonds. Then,
    for every bond that is a triple bond connecting a carbon and a nitrogen, 
    it checks if the carbon has exactly two neighbors (once H's are explicit) and if 
    its non-nitrogen neighbor is not a hydrogen or a metal. Only such groups are accepted.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains at least one valid substituted nitrile group.
        str: An explanation for the classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that we can properly count the bonds at nitrile carbons.
    mol = Chem.AddHs(mol)
    
    valid_nitrile_found = False
    # Iterate over all bonds in the molecule.
    for bond in mol.GetBonds():
        # Look for triple bonds.
        if bond.GetBondType() == Chem.BondType.TRIPLE:
            # Identify the two atoms.
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # We want one atom to be carbon and the other to be nitrogen.
            if (atom1.GetSymbol() == "C" and atom2.GetSymbol() == "N"):
                nitrile_c = atom1
                nitrile_n = atom2
            elif (atom2.GetSymbol() == "C" and atom1.GetSymbol() == "N"):
                nitrile_c = atom2
                nitrile_n = atom1
            else:
                continue  # not a C#N group
            
            # For a proper nitrile the C should be sp-hybridized with exactly 2 neighbors.
            # (One neighbor is the nitrile N, and the other should be the substituent.)
            if nitrile_c.GetDegree() != 2:
                continue
            
            # Get all neighbors of the nitrile carbon.
            neighbors = nitrile_c.GetNeighbors()
            # Identify the neighbor that is not the nitrile nitrogen.
            substituent = None
            for nb in neighbors:
                if nb.GetIdx() == nitrile_n.GetIdx():
                    continue
                substituent = nb
                break
                
            # If for some reason there is no substituent, skip.
            if substituent is None:
                continue
            
            # If the substituent atom is hydrogen, then this is just HC#N.
            if substituent.GetAtomicNum() == 1:
                continue
            
            # If the substituent atom is a metal, then reject this nitrile group.
            if substituent.GetSymbol() in metal_symbols:
                continue
            
            # If we reach here, we have a nitrile carbon (C#N) where the nitrile carbon
            # has a substituent (other than a hydrogen or metal). This is classified as RC#N.
            valid_nitrile_found = True
            return True, "Molecule contains a substituted nitrile group (RC#N)"
    
    # If no valid nitrile groups were found.
    return False, "No substituted nitrile group (RC#N) found; only HCN-like or metal‐bound nitrile(s) detected"

# Example usage:
if __name__ == "__main__":
    # Test some SMILES strings:
    test_smiles_list = [
        "CC#N",                    # acetonitrile (valid nitrile)
        "N#C[Fe]C#N",              # iron dicyanide (should not be counted)
        "C#N",                     # hydrogen cyanide (invalid as not substituted)
    ]
    for smi in test_smiles_list:
        result, reason = is_nitrile(smi)
        print(f"{smi}: {result} -- {reason}")