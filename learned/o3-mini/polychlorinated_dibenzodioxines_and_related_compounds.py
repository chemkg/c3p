"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Polychlorinated dibenzodioxines and related compounds
These are organochlorine (and/or organobromine) compounds that include
polychlorinated dibenzodioxines, polychlorinated dibenzofurans, and
polychlorinated/polybrominated biphenyls. Such compounds typically have a
core aromatic structure (dioxin, dibenzofuran, or biphenyl) and multiple
halogen substituents.
"""

from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to polychlorinated dibenzodioxines and
    related compounds based on its SMILES string. This includes checking for
    aromatic cores (like dibenzodioxin, dibenzofuran or biphenyl motifs) and the
    presence of several chlorine and/or bromine atoms.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule belongs to this class, False otherwise
        str: Reason for classification decision
    """
    # Parse SMILES string into molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule has computed aromaticity.
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
    
    # Count halogen atoms (chlorine and bromine).
    halogen_atoms = {17, 35}  # atomic numbers for Cl and Br respectively.
    halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in halogen_atoms)
    if halogen_count < 2:
        return False, "Not enough halogen atoms (Cl/Br) for polyhalogenation"

    # Check that the molecule has at least 2 aromatic rings.
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        # Verify every atom in the ring is aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
    if aromatic_ring_count < 2:
        return False, "Not enough aromatic rings to support a dioxin/dibenzofuran/biphenyl structure"
    
    # Define SMARTS for typical aromatic cores.
    # 1. Dioxin core: two benzene rings connected by two oxygen atoms.
    dioxin_smarts = "O1c2ccccc2Oc3ccccc31"
    dioxin_pattern = Chem.MolFromSmarts(dioxin_smarts)
    
    # 2. Dibenzofuran core: two benzene rings connected by one oxygen atom.
    dibenzofuran_smarts = "c1ccc2occc2c1"
    dibenzofuran_pattern = Chem.MolFromSmarts(dibenzofuran_smarts)
    
    # 3. Biphenyl core: two benzene rings directly connected.
    biphenyl_smarts = "c1ccccc1-c2ccccc2"
    biphenyl_pattern = Chem.MolFromSmarts(biphenyl_smarts)
    
    # Check for at least one of the cores using substructure search.
    has_dioxin = mol.HasSubstructMatch(dioxin_pattern)
    has_dibenzofuran = mol.HasSubstructMatch(dibenzofuran_pattern)
    has_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)
    
    if not (has_dioxin or has_dibenzofuran or has_biphenyl):
        return False, "Does not contain a dioxin, dibenzofuran, or biphenyl core"
    
    # Optionally, we can be more strict by requiring more than a minimal number of halogen atoms.
    # Here we use at least 3 halogen atoms, for example.
    if halogen_count < 3:
        return False, "Insufficient halogen substitutions for a polychlorinated structure"
    
    # If all conditions are met, the molecule classifies as polychlorinated
    # dibenzodioxines and related compounds.
    return True, "Contains aromatic core with multiple halogen atoms consistent with the target class"

# (Optional) For testing purposes:
if __name__ == "__main__":
    # Test examples: uncomment to try
    examples = [
        ("1,2,3,7,8-Pentachlorodibenzodioxin", "Clc1cc2Oc3cc(Cl)c(Cl)c(Cl)c3Oc2cc1Cl"),
        ("2,2',3,3',5,5'-hexachlorobiphenyl", "Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)cc(Cl)c1Cl"),
        ("3,5-dichlorobiphenyl", "Clc1cc(Cl)cc(c1)-c1ccccc1")
    ]
    for name, smi in examples:
        result, reason = is_polychlorinated_dibenzodioxines_and_related_compounds(smi)
        print(f"{name}: {result} ({reason})")