"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: Organohalogen compound
Definition: A compound containing at least one carbon-halogen bond (halogen = F, Cl, Br, I).
"""

from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    Instead of using a single SMARTS query, this function first isolates the main (largest)
    fragment of the molecule and then iterates over all bonds to check for a bond
    between a carbon and a halogen (F, Cl, Br, I).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule (main fragment) contains at least one carbon-halogen bond,
              False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # In cases where the molecule contains multiple disconnected fragments
    # (e.g., an organic structure plus inorganic counterions), use the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True)
    main_mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # Define the set of halogen atomic numbers: F (9), Cl (17), Br (35), I (53)
    halogen_atomic_numbers = {9, 17, 35, 53}
    
    # Iterate over all bonds in the main fragment.
    # Check each bond: if one atom is carbon (atomic number 6) and the other is a halogen.
    for bond in main_mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() in halogen_atomic_numbers) or \
           (atom2.GetAtomicNum() == 6 and atom1.GetAtomicNum() in halogen_atomic_numbers):
            return True, "Molecule contains a carbon-halogen bond on its main organic fragment"
    
    return False, "No carbon-halogen bond found in the main organic fragment"

# Example usage for testing (can be commented out in production)
if __name__ == "__main__":
    test_smiles = [
        "COc1ccc(cc1)C(=O)C(Br)CS(=O)(=O)c1ccc(C)cc1",  # 2-bromo-1-(4-methoxyphenyl)-3-[(4-methylphenyl)sulfonyl]-1-propanone
        "Oc1c(Cl)cc(Cl)cc1-c1ccccc1",                     # 2-Hydroxy-3,5-dichlorobiphenyl
        "C(O)(C(=O)CBr)=O",                              # 3-bromopyruvic acid
        "CCO",                                          # Ethanol, no halogen
    ]
    
    for smi in test_smiles:
        result, reason = is_organohalogen_compound(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")