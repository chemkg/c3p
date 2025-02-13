"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: Organohalogen compound
Definition: A compound containing at least one carbon-halogen bond (halogen = F, Cl, Br, I).
Improved approach: Select the main organic fragment as the one with the highest number of carbon atoms.
"""

from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule (main organic fragment) is an organohalogen compound based on its SMILES string.
    The function first filters out fragments lacking carbon atoms, then selects the fragment with the highest carbon count.
    Finally, it iterates over all bonds in that fragment to check for a bond between a carbon atom and a halogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the main organic fragment contains at least one carbon-halogen bond, False otherwise.
        str: Explanation for the classification.
    """
    # Try to parse the SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all disconnected fragments as separate molecules.
    fragments = Chem.GetMolFrags(mol, asMols=True)
    
    # Filter fragments to include only those with at least one carbon atom.
    organic_frags = []
    for frag in fragments:
        if any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms()):
            organic_frags.append(frag)
    if not organic_frags:
        return False, "No organic fragment (with carbon atoms) found"
    
    # Select the fragment with the highest count of carbon atoms.
    def carbon_count(m):
        return sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 6)
    main_mol = max(organic_frags, key=carbon_count)
    
    # Define the halogen atomic numbers: F (9), Cl (17), Br (35), I (53)
    halogen_atomic_nums = {9, 17, 35, 53}
    
    # Iterate over all bonds in the main fragment.
    for bond in main_mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() in halogen_atomic_nums) or \
           (atom2.GetAtomicNum() == 6 and atom1.GetAtomicNum() in halogen_atomic_nums):
            return True, "Molecule contains a carbon-halogen bond on its main organic fragment"
    
    return False, "No carbon-halogen bond found in the main organic fragment"

# Example usage for testing (can be commented out in production):
if __name__ == "__main__":
    test_smiles = [
        "COc1ccc(cc1)C(=O)C(Br)CS(=O)(=O)c1ccc(C)cc1",    # 2-bromo-1-(4-methoxyphenyl)-3-[(4-methylphenyl)sulfonyl]-1-propanone
        "Oc1c(Cl)cc(Cl)cc1-c1ccccc1",                       # 2-Hydroxy-3,5-dichlorobiphenyl
        "C(O)(C(=O)CBr)=O",                                # 3-bromopyruvic acid
        "C1=CC=C2C(=C1)C(=O)C3=CC=CC=C3C2(CC4=CC(=NC=C4)F)CC5=CC(=NC=C5)F",  # a false positive challenge case
        "Brn1ccc2ccccc12",                                # 1-bromoindole; challenge case
    ]
    
    for smi in test_smiles:
        result, reason = is_organohalogen_compound(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")