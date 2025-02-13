"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: Organoiodine compounds
Definition: An organoiodine compound is defined here as an organic molecule (i.e. no metals or other unexpected elements)
that contains at least one “normal” carbon–iodine bond – namely a single (or aromatic) bond between a neutral C and iodine –
and whose overall iodine content is not excessively high relative to the number of carbons.
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    For our purposes, an organoiodine compound is defined as an organic molecule (i.e. containing only expected non‐metal atoms)
    that has at least one valid (neutral) carbon–iodine bond and does not show an unusually high iodine density.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an organoiodine compound, False otherwise.
        str: A reason string for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, check that the molecule does not contain any atoms outside a typical organic set.
    # Here we allow H, B, C, N, O, F, Si, P, S, Cl, Br, and I.
    allowed_atoms = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in allowed_atoms:
            return False, f"Molecule contains disallowed atom with atomic number {atomic_num} (likely a metal or inorganic element)"
    
    # Count carbons and iodines in the molecule – later we can check their ratio.
    carbon_count = 0
    iodine_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
        elif atom.GetAtomicNum() == 53:
            iodine_count += 1

    if iodine_count == 0:
        return False, "No iodine atoms present in the molecule"
    if carbon_count == 0:
        return False, "No carbon atoms present in the molecule"
    
    # If the iodine content is unusually high relative to carbons, this may indicate a salt/contrast agent.
    ratio = iodine_count / carbon_count
    # Here we choose a cutoff (e.g. 0.30). (This value can be tuned; many genuine organoiodine compounds have a lower ratio.)
    if ratio > 0.30:
        return False, f"Iodine/carbon ratio too high ({iodine_count}/{carbon_count} = {ratio:.2f}); not a typical organoiodine compound"
    
    # Now, look for a valid carbon–iodine bond.
    # We iterate over all bonds in the molecule and search for one where one atom is C and the other is I.
    valid_ci_found = False
    for bond in mol.GetBonds():
        # We accept bonds that are either SINGLE or aromatic.
        bond_type = bond.GetBondType()
        if bond_type not in (Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.AROMATIC):
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Check one atom is carbon (6) and the other is iodine (53); also require both atoms to have formal charge 0.
        if ((a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 53) or (a1.GetAtomicNum() == 53 and a2.GetAtomicNum() == 6)):
            if a1.GetFormalCharge() == 0 and a2.GetFormalCharge() == 0:
                valid_ci_found = True
                break

    if valid_ci_found:
        return True, "Molecule contains at least one valid carbon–iodine bond in an organic context"
    else:
        return False, "No valid carbon–iodine bond found in the molecule"

# Example test calls (you can remove or comment these out in actual module usage)
if __name__ == '__main__':
    # True positive example: 3-(2-iodacetamido)-PROXYL
    smiles_tp = "N1([O])C(C(CC1(C)C)NC(CI)=O)(C)C"
    is_org, reason = is_organoiodine_compound(smiles_tp)
    print(smiles_tp, is_org, reason)
    
    # False positive example: Calcium diiodostearate (should return False because of metal content)
    smiles_fp = "[Ca+2].C(CCCCCCCC(CCC(CCCCCC)I)I)([O-])=O.C(CCCCCCCC(CCC(CCCCCC)I)I)([O-])=O"
    res, reason_fp = is_organoiodine_compound(smiles_fp)
    print(smiles_fp, res, reason_fp)