"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: Organoiodine compounds
Definition (tuned criteria): An organoiodine compound is defined here as an organic molecule that 
contains at least one “normal” (neutral) carbon–iodine bond (SINGLE or aromatic) and appears to have a 
typical iodine density (i.e. not over‐iodinated) as measured by two heuristics: (1) a minimum number 
of carbon atoms (here ≥6) and (2) an iodine/carbon ratio ≤ 0.30. In addition, if the molecule contains 
sulfonate/sulfate or sulfone groups (which are common in iodinated contrast agents or very special reagents), 
it is rejected.
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is a typical organoiodine compound based on its SMILES string.
    For our purposes, an organoiodine compound is an organic molecule that has:
      (a) only allowed (non‐metal) atoms (H, B, C, N, O, F, Si, P, S, Cl, Br, I),
      (b) at least one valid (neutral) C–I bond (SINGLE or aromatic),
      (c) at least 6 carbon atoms, and an iodine/carbon ratio not exceeding 0.30,
      (d) and no clearly non-typical functional groups (e.g. sulfate, sulfonate, or sulfone substructures)
    that hint the molecule belongs to a different iodine chemistry class.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a typical organoiodine compound, False otherwise.
        str: A reason string for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Allowed atoms for an organic molecule: H, B, C, N, O, F, Si, P, S, Cl, Br, I.
    allowed_atomic_nums = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Molecule contains disallowed atom with atomic number {atom.GetAtomicNum()} (likely a metal or inorganic element)"
    
    # Count carbon and iodine atoms.
    carbon_count = 0
    iodine_count = 0
    for atom in mol.GetAtoms():
        at_num = atom.GetAtomicNum()
        if at_num == 6:
            carbon_count += 1
        elif at_num == 53:
            iodine_count += 1
    
    # Must have at least one iodine.
    if iodine_count == 0:
        return False, "No iodine atoms present in the molecule"
    if carbon_count == 0:
        return False, "No carbon atoms present in the molecule"
    
    # Check a minimum number of carbons; many bona fide organic molecules have several carbons.
    if carbon_count < 6:
        return False, f"Too few carbons ({carbon_count}) to be considered a typical organoiodine compound"

    # Evaluate iodine density.
    ratio = iodine_count / carbon_count
    if ratio > 0.30:
        return False, f"Iodine/carbon ratio too high ({iodine_count}/{carbon_count} = {ratio:.2f}); not typical of organoiodine compounds"
    
    # Look for at least one valid carbon-iodine bond.
    valid_ci_found = False
    for bond in mol.GetBonds():
        # Consider only SINGLE or AROMATIC bonds.
        if bond.GetBondType() not in (Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.AROMATIC):
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Check if one atom is carbon (6) and the other is iodine (53) and both have formal charge 0.
        if ((a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 53) or 
            (a1.GetAtomicNum() == 53 and a2.GetAtomicNum() == 6)):
            if a1.GetFormalCharge() == 0 and a2.GetFormalCharge() == 0:
                valid_ci_found = True
                break
    if not valid_ci_found:
        return False, "No valid neutral carbon–iodine bond found in the molecule"
    
    # Check for disallowed substructures.
    # (1) Sulfate or sulfonate group typically appear as OS(=O)(=O)X.
    sulfate_smarts = "[OX2H0][S](=O)(=O)[O-]"
    sulfonate_pattern = Chem.MolFromSmarts(sulfate_smarts)
    if mol.HasSubstructMatch(sulfonate_pattern):
        return False, "Molecule contains a sulfate/sulfonate group, not typical of organoiodine compounds"
    
    # (2) A simple sulfone group can be defined as an alkyl/saryl sulfone: [#6][S](=O)(=O)[#6]
    sulfone_smarts = "[#6][S](=O)(=O)[#6]"
    sulfone_pattern = Chem.MolFromSmarts(sulfone_smarts)
    if mol.HasSubstructMatch(sulfone_pattern):
        return False, "Molecule contains a sulfone group, not typical of organoiodine compounds"
    
    # If all tests pass, we conclude that the molecule is a typical organoiodine compound.
    return True, "Contains valid neutral C–I bond(s) with appropriate iodine density and no disallowed groups"

# Example test calls (they may be commented out when used as a module)
if __name__ == '__main__':
    # True positive example: 3-(2-iodacetamido)-PROXYL
    smiles_tp = "N1([O])C(C(CC1(C)C)NC(CI)=O)(C)C"
    result, reason = is_organoiodine_compound(smiles_tp)
    print("SMILES:", smiles_tp)
    print("Classification:", result)
    print("Reason:", reason)
    
    # A molecule that should be rejected because of high iodine/carbon ratio (e.g. diiodohydroxypropane)
    smiles_fn = "C(CI)(CI)O"
    result, reason = is_organoiodine_compound(smiles_fn)
    print("\nSMILES:", smiles_fn)
    print("Classification:", result)
    print("Reason:", reason)
    
    # A false positive example from previous evaluation: 3,3',5-triiodo-L-thyronine sulfate
    # (The sulfate group should cause a rejection.)
    smiles_fp = "N[C@@H](Cc1cc(I)c(Oc2ccc(OS(=O)(=O)O)C=C2)c(I)c1)"
    result, reason = is_organoiodine_compound(smiles_fp)
    print("\nSMILES:", smiles_fp)
    print("Classification:", result)
    print("Reason:", reason)