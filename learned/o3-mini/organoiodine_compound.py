"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: Organoiodine compounds
Definition: An organoiodine compound is defined as an organic molecule that contains at least one neutral carbon–iodine bond.
Additional filters (to reject clearly non‐organic or contrast-agent compounds) include:
  - Only allowed atoms are used (H, B, C, N, O, F, Si, P, S, Cl, Br, I).
  - Molecules containing sulfate/sulfonate or sulfone substructures are rejected.
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound (i.e. containing at least one neutral carbon–iodine bond)
    while rejecting molecules with disallowed atoms and substructures (e.g. sulfate/sulfonate or sulfone groups).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an organoiodine compound, False otherwise.
        str: A reason message for the classification decision.
    """
    # Convert SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define allowed atomic numbers for a typical organic molecule.
    allowed_atomic_nums = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}  # H, B, C, N, O, F, Si, P, S, Cl, Br, I
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Molecule contains disallowed atom with atomic number {atom.GetAtomicNum()} (likely a metal or inorganic element)"

    # Check for disallowed substructures.
    # (1) Sulfate or sulfonate group: e.g. OS(=O)(=O)[O-]
    sulfate_smarts = "[OX2][S](=O)(=O)[O-]"
    sulfate_pattern = Chem.MolFromSmarts(sulfate_smarts)
    if mol.HasSubstructMatch(sulfate_pattern):
        return False, "Molecule contains a sulfate/sulfonate group, not typical of the desired organoiodine class"
    
    # (2) Simple sulfone group: e.g. [#6][S](=O)(=O)[#6]
    sulfone_smarts = "[#6][S](=O)(=O)[#6]"
    sulfone_pattern = Chem.MolFromSmarts(sulfone_smarts)
    if mol.HasSubstructMatch(sulfone_pattern):
        return False, "Molecule contains a sulfone group, not typical of the desired organoiodine class"

    # Loop over bonds to find at least one valid (neutral) carbon–iodine bond.
    # We restrict our search to SINGLE and AROMATIC bonds.
    for bond in mol.GetBonds():
        if bond.GetBondType() not in (Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.AROMATIC):
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Check if one atom is carbon (6) and the other is iodine (53)
        if ((a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 53) or 
            (a1.GetAtomicNum() == 53 and a2.GetAtomicNum() == 6)):
            # Both atoms should be neutral.
            if a1.GetFormalCharge() == 0 and a2.GetFormalCharge() == 0:
                return True, "Contains at least one neutral carbon–iodine bond"
                
    return False, "No valid neutral carbon–iodine bond found in the molecule"

# Example testing when running as a standalone script
if __name__ == '__main__':
    # Example: 3-(2-iodacetamido)-PROXYL; should be classified as an organoiodine compound.
    smiles_tp = "N1([O])C(C(CC1(C)C)NC(CI)=O)(C)C"
    result, reason = is_organoiodine_compound(smiles_tp)
    print("SMILES:", smiles_tp)
    print("Classification:", result)
    print("Reason:", reason)
    
    # Example: diiodohydroxypropane; although small, it contains two C–I bonds.
    smiles_small = "C(CI)(CI)O"
    result, reason = is_organoiodine_compound(smiles_small)
    print("\nSMILES:", smiles_small)
    print("Classification:", result)
    print("Reason:", reason)
    
    # Example: 3,3',5-triiodo-L-thyronine sulfate; this contains a sulfate group and should be rejected.
    smiles_fp = "N[C@@H](Cc1cc(I)c(Oc2ccc(OS(=O)(=O)O)C=C2)c(I)c1)"
    result, reason = is_organoiodine_compound(smiles_fp)
    print("\nSMILES:", smiles_fp)
    print("Classification:", result)
    print("Reason:", reason)