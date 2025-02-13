"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: Primary Amine
Definition: A compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.
For this classifier, the molecule is considered a primary amine if it contains at least one primary amine group.
A valid primary amine group is defined as an –NH2 group (with exactly two hydrogens) bonded to one non-hydrogen atom 
(that is, derived from ammonia by replacing one hydrogen by a hydrocarbyl group) and not being involved in an amide bond.
Also, if the molecule contains any amide bonds (indicative of peptides or similar compounds), the molecule is rejected.
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine (R–NH2) is derived from ammonia by replacing one hydrogen with a hydrocarbyl group.
    For our purposes, the molecule is considered a primary amine if it meets the following criteria:
      - It does NOT contain any amide bonds (the presence of a pattern "C(=O)N" is grounds for rejection).
      - It contains at least one primary amine group. 
        A primary amine group is found by matching the SMARTS pattern "[NX3;H2;!$(NC=O)]":
            • NX3: a trivalent nitrogen,
            • H2: exactly two hydrogens attached,
            • !$(NC=O): not directly bonded to a carbonyl carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a primary amine, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are accessible.
    mol = Chem.AddHs(mol)
    
    # First, reject molecules that contain an amide bond.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    if amide_smarts is None:
        return False, "Error generating amide SMARTS pattern"
    if mol.HasSubstructMatch(amide_smarts):
        return False, "Molecule contains amide bond(s) typical of peptides/amides"
    
    # Define the SMARTS for a primary amine:
    # [NX3;H2] means trigonal nitrogen with exactly two hydrogens.
    # !$(NC=O) ensures that the nitrogen is not bonded to a carbonyl group.
    primary_amine_smarts = Chem.MolFromSmarts("[NX3;H2;!$(NC=O)]")
    if primary_amine_smarts is None:
        return False, "Error generating primary amine SMARTS pattern"
    
    # Look for a match of the primary amine group in the molecule.
    if mol.HasSubstructMatch(primary_amine_smarts):
        return True, "Molecule contains a valid primary amine (R–NH2) group."
    
    return False, "No valid primary amine group (with two hydrogens and one substituent not involved in an amide) found."

# Example test calls (you may remove these before deployment)
if __name__ == "__main__":
    test_smiles_list = [
        # acid fuchsin (free acid form)
        "Cc1cc(cc(c1N)S(O)(=O)=O)C(=C1\\C=CC(=N)C(=C1)S(O)(=O)=O)\\c1ccc(N)c(c1)S(O)(=O)=O",
        # N-[(3-methoxyphenyl)methyl]-N-methyl-1-[1-[2-(2-methylphenyl)ethyl]-3-piperidinyl]methanamine
        "CC1=CC=CC=C1CCN2CCCC(C2)CN(C)CC3=CC(=CC=C3)OC",
        # (R)-clenbuterol
        "CC(C)(C)NC[C@H](O)c1cc(Cl)c(N)c(Cl)c1",
        # 2-Methylbutylamine
        "NC[C@H](CC)C",
        # clenbuterol
        "CC(C)(C)NCC(O)c1cc(Cl)c(N)c(Cl)c1",
        # Phenelzine
        "NNCCc1ccccc1",
        # methylamine (a simple primary amine)
        "CN",
        # aniline
        "Nc1ccccc1",
        # pentan-1-amine
        "CCCCCN"
    ]
    for smi in test_smiles_list:
        result, reason = is_primary_amine(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")