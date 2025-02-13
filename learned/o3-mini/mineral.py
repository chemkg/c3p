"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: minerals
A mineral is normally an inorganic crystalline (or formerly crystalline) substance,
often formed by geological processes. Heuristics used:
  - Many minerals are composed of ionic fragments (multiple disconnected components).
  - Many minerals include metal atoms, but if the molecule contains many carbons and hydrogens,
    it likely is an organic salt.
  - Many show formal charges.
  - Many inorganic minerals either lack carbon entirely or have very few C atoms.
  
Note: Because crystallinity or the geological history is not encoded in SMILES,
this classifier uses structural proxies and may not cover all edge cases.
"""

from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is likely a mineral based on its SMILES string using
    several heuristic checks.
    
    Heuristics:
      - A fragment count >1 (disconnected parts) is typical of ionic compounds.
      - Presence of metal atoms in a molecule that is not clearly organic.
      - Molecules with formal charges that are not clearly organic.
      - Molecules containing no carbon (often seen in inorganic minerals).
      
    To better differentiate organic compounds from inorganic salts, we add explicit hydrogens
    and consider a molecule as "organic" if it contains at least 3 carbon atoms and 3 hydrogen atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a mineral, False otherwise.
        str: Explanation for the decision.
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to get a better count of hydrogen atoms (RDKit stores H implicitly)
    mol_with_H = Chem.AddHs(mol)
    
    # Calculate the number of fragments. Many minerals (salts) come as disconnected ions.
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
    
    # Initialize counters and flags for the heuristics
    metal_found = False
    has_formal_charge = False
    count_C = 0  # count of carbon atoms (using original molecule which has all atoms)
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol == "C":
            count_C += 1
        if symbol in metal_symbols:
            metal_found = True
        # Check if any atom has a non-zero formal charge
        if atom.GetFormalCharge() != 0:
            has_formal_charge = True
            
    # Count explicit hydrogens using the molecule with hydrogens
    count_H = sum(1 for atom in mol_with_H.GetAtoms() if atom.GetSymbol() == "H")
    
    # Define a heuristic to decide if a molecule is "organic"
    # Here, we assume a molecule is organic if it contains at least 3 carbons and 3 hydrogens.
    is_organic = (count_C >= 3 and count_H >= 3)
    
    # Criterion 1: If the molecule is composed of multiple disconnected fragments,
    # we assume an ionic (and thus potentially mineral) nature.
    if multiple_fragments:
        return True, "Contains multiple disconnected fragments, consistent with an ionic compound typical of minerals."
    
    # Criterion 2: If the molecule contains no carbon,
    # it is likely a purely inorganic substance.
    if count_C == 0:
        return True, "Molecule contains no carbon, a hallmark of many inorganic minerals."
    
    # Criterion 3: If the molecule has metal(s) and is not clearly organic,
    # then classify it as a mineral.
    if metal_found and not is_organic:
        return True, "Contains metal element(s) typical of minerals and is not obviously organic."
    
    # Criterion 4: If the molecule shows formal charges and is not clearly organic,
    # classify as mineral.
    if has_formal_charge and not is_organic:
        return True, "Molecule contains formal charges, suggesting an ionic nature typical of minerals, and is not obviously organic."
    
    # If none of the heuristics for minerality are met,
    # then we treat the molecule as not a mineral.
    return False, "Does not exhibit typical inorganic mineral characteristics (e.g. disconnected ionic fragments, low carbon count, or metal ions in a non-organic framework)."

# Example usage (for testing):
if __name__ == '__main__':
    test_smiles_list = [
        "[Na+].[Na+].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]OB1O[B-]2(O[H])OB(O[H])O[B-](O[H])(O1)O2",  # disodium tetraborate octahydrate
        "[Fe++].[S-][S-]",  # pyrite
        "[OH-].[OH-].[Ca++]",  # calcium dihydroxide
        "[Na+].CCCC([O-])=O",  # sodium butyrate (intended to be organic)
    ]
    for smi in test_smiles_list:
        result, reason = is_mineral(smi)
        print("SMILES:", smi)
        print("is_mineral:", result)
        print("Reason:", reason)
        print("------")