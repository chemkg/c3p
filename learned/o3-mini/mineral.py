"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: minerals
A mineral is normally an inorganic crystalline (or formerly crystalline) substance,
often formed by geological processes. Because properties such as crystallinity or geological history
are not encoded in a SMILES string, this classifier uses structural proxies.
Heuristics used include:
  - Many minerals are ionic compounds that appear as disconnected fragments.
  - Many minerals contain metal atoms or formal charges.
  - Many inorganic minerals either have no carbon or only non-organic carbon skeletons.
We now also check individual fragments: if any fragment has at least 3 carbons and 3 hydrogens (after adding H’s),
it is taken to be organic, and the overall molecule is not classified as a mineral.
Additionally, very small species (fewer than three atoms) are not classified as minerals.
"""

from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is likely a mineral based on its SMILES string using
    several heuristic checks.
    
    Heuristics used:
      - If the molecule is very small (fewer than 3 atoms), do not call it a mineral.
      - Disconnected fragments are examined individually. If any fragment qualifies as organic 
        (defined as having at least 3 carbon atoms and 3 hydrogen atoms, considering implicit hydrogens),
        then the molecule is classified as non-mineral (an organic salt) despite ionic disconnection.
      - If no organic fragment is present then:
             * A molecule with more than one fragment is classified as a mineral.
             * In a single fragment, if there is no carbon or if metal(s) or formal charges are present,
               we classify it as mineral.
    
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
    
    # Reject very tiny species: a realistic mineral should have at least 3 atoms.
    if mol.GetNumAtoms() < 3:
        return False, "Molecule too small to be considered a mineral"
    
    # Define a set of metal element symbols that are common in minerals.
    metal_symbols = {
        "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
        "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
        "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce",
        "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
        "Po", "Fr", "Ra", "Ac", "Th", "U", "Np", "Pu"
    }
    
    # Helper function to check if a fragment is clearly organic.
    def is_fragment_organic(frag):
        # Add explicit hydrogens to count H atoms reliably.
        frag_with_H = Chem.AddHs(frag)
        count_C = sum(1 for atom in frag_with_H.GetAtoms() if atom.GetSymbol() == "C")
        count_H = sum(1 for atom in frag_with_H.GetAtoms() if atom.GetSymbol() == "H")
        # Define a fragment as organic if it has at least 3 carbons and at least 3 hydrogens.
        return (count_C >= 3 and count_H >= 3)
    
    # Get individual fragments (ions/molecular pieces) from the molecule.
    fragments = Chem.GetMolFrags(mol, asMols=True)
    
    # Examine if any fragment is clearly organic.
    for frag in fragments:
        if is_fragment_organic(frag):
            return False, ("Contains organic fragment(s) (at least 3 C and 3 H atoms in a fragment), "
                           "indicating an organic salt rather than a typical inorganic mineral.")
    
    # If no organic fragments are present, then use additional heuristics.
    # Criterion 1: If there are multiple disconnected fragments, the ionic nature is consistent with minerals.
    if len(fragments) > 1:
        return True, "Contains multiple disconnected fragments, consistent with an ionic compound typical of minerals."
    
    # For a single-fragment molecule, further examine its properties.
    count_C = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")
    # Criterion 2: Molecules with no carbon are often inorganic minerals.
    if count_C == 0:
        return True, "Molecule contains no carbon, a hallmark of many inorganic minerals."
    
    # Check for metals and any formal charges in the molecule.
    metal_found = any(atom.GetSymbol() in metal_symbols for atom in mol.GetAtoms())
    has_formal_charge = any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms())
    
    # Criterion 3: If the molecule contains metal(s) or formal charges (and is not organic by the previous test),
    # then classify it as a mineral.
    if metal_found or has_formal_charge:
        return True, "Contains metal element(s) or formal charges typical of inorganic minerals."
    
    # If none of the inorganic heuristics match, then it is not classified as a mineral.
    return False, "Does not exhibit typical inorganic mineral characteristics."

# Example usage (for testing):
if __name__ == '__main__':
    test_smiles_list = [
        "[Na+].[Na+].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]OB1O[B-]2(O[H])OB(O[H])O[B-](O[H])(O1)O2",  # disodium tetraborate octahydrate
        "[Fe++].[S-][S-]",  # pyrite
        "[OH-].[OH-].[Ca++]",  # calcium dihydroxide
        "[O-]S([O-])(=O)=O.[Ba+2]",  # barium sulfate
        "[Na+].CCCC([O-])=O",  # sodium butyrate (should be classified as non-mineral)
        "[7Li]",  # lithium-7 atom (too small to be a mineral)
    ]
    for smi in test_smiles_list:
        result, reason = is_mineral(smi)
        print("SMILES:", smi)
        print("is_mineral:", result)
        print("Reason:", reason)
        print("------")