"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: Organohalogen compound
Definition: A compound is considered an organohalogen compound if at least one major organic fragment (one in which the carbon count is at least 50% of the largest fragment’s carbon count) contains a carbon-halogen bond (with halogen = F, Cl, Br, I).
Note: This approach first splits the molecule into fragments containing at least one carbon; then it selects those fragments that are “major” (by carbon count) and finally checks for a bond where one atom is carbon and the other one is F, Cl, Br, or I.
"""

from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound.
    The function splits the molecule into fragments that contain carbon atoms.
    It then identifies the major fragments – those with a carbon count at least 50% of the maximum carbon count among fragments.
    Finally, it checks those fragments for the presence of a carbon-halogen bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if at least one major organic fragment contains a C–X bond (X = F, Cl, Br, I), False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Split into disconnected fragments.
    frags = Chem.GetMolFrags(mol, asMols=True)
    # Filter fragments containing at least one carbon atom.
    organic_frags = []
    for f in frags:
        if any(atom.GetAtomicNum() == 6 for atom in f.GetAtoms()):
            organic_frags.append(f)
    if not organic_frags:
        return False, "No organic fragment (containing carbon) found"
    
    # Count carbon atoms in each fragment.
    def carbon_count(m):
        return sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 6)
    frag_carbon_counts = [carbon_count(f) for f in organic_frags]
    max_carbons = max(frag_carbon_counts)
    
    # Decide which fragments are "major" (at least 50% as many carbons as the largest)
    major_frags = [frag for frag, cnt in zip(organic_frags, frag_carbon_counts) if cnt >= 0.5 * max_carbons]
    if not major_frags:
        return False, "No major organic fragment found"
    
    # Define halogen atomic numbers.
    halogen_atomic_nums = {9, 17, 35, 53}
    
    # Check each bond in each major fragment.
    for frag in major_frags:
        for bond in frag.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Only count a bond as a C–X bond if one atom is carbon and the other is one of the defined halogens.
            if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() in halogen_atomic_nums) or \
               (a2.GetAtomicNum() == 6 and a1.GetAtomicNum() in halogen_atomic_nums):
                return True, "Major organic fragment contains a carbon-halogen bond"
    
    return False, "No carbon-halogen bond found in the major organic fragment"

# Example usage (for testing; remove or comment out when integrating into production):
if __name__ == "__main__":
    test_smiles_list = [
        "COc1ccc(cc1)C(=O)C(Br)CS(=O)(=O)c1ccc(C)cc1",  # True positive
        "C=1(C=C(N(N1)CC(=O)N2CCC(CC2)C3=NC(=CS3)C=4C[C@@](ON4)(C=5C(=CC=CC5Cl)OS(C)(=O)=O)[H])C(F)F)C(F)F",  # True positive
        "Oc1c(Cl)cc(Cl)cc1-c1ccccc1",  # True positive
        "C(O)(C(=O)CBr)=O",  # True positive
        "C1=CC=C2C(=C1)C(=O)C3=CC=CC=C3C2(CC4=CC(=NC=C4)F)CC5=CC(=NC=C5)F",  # Challenge false positive case
        "Brn1ccc2ccccc12",  # 1-bromoindole; expected false because Br is not bound to carbon
        "Cln1c(=O)n(Cl)c(=O)n(Cl)c1=O",  # trichloro-triazinone; expected false (no carbon)
    ]
    
    for smi in test_smiles_list:
        result, reason = is_organohalogen_compound(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")