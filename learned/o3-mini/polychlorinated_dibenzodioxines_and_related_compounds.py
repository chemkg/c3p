"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Polychlorinated dibenzodioxines and related compounds.
These include organochlorine/organobromine compounds such as 
polychlorinated dibenzodioxins, dibenzofurans, and polyhalogenated biphenyls.
The approach is to identify a core aromatic motif (dioxin, dibenzofuran or biphenyl)
and then count the number of halogen substituents (Cl or Br) attached to that core.
A minimum of 2 such substituents is required.
"""

from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to polychlorinated dibenzodioxines and related compounds.
    
    The function first sanitizes the molecule, then looks for one of the 
    following aromatic cores by substructure matching:
      - Dibenzodioxin core (two benzene rings connected by two oxygen bridges)
      - Dibenzofuran core (two benzene rings connected by one oxygen)
      - Biphenyl core (two benzene rings directly connected)
      
    For every match of one of these cores, only the halogen atoms that are directly
    attached (neighbors but not part of the core itself) are counted.
    
    If at least one core is found that has >=2 attached Cl or Br atoms, the molecule
    is classified as being in this chemical class.
    
    If none of the cores meet the criteria, the function returns False with an explanation.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule belongs to this class, otherwise False.
        str: Explanation for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    try:
        # Sanitize the molecule (this ensures aromaticity is computed etc).
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"
    
    # Define atomic numbers for considered halogens.
    halogen_znums = {17, 35}  # Cl and Br

    # Define SMARTS patterns for the target aromatic cores.
    # 1. Dibenzodioxin (dioxin) core: two benzene rings connected by two oxygen bridges.
    dioxin_smarts = "c1ccc2c(c1)Oc3ccccc3O2"
    dioxin_core = Chem.MolFromSmarts(dioxin_smarts)
    
    # 2. Dibenzofuran core: two benzene rings connected by one oxygen.
    dibenzofuran_smarts = "c1ccc2occc2c1"
    dibenzofuran_core = Chem.MolFromSmarts(dibenzofuran_smarts)
    
    # 3. Biphenyl core: two benzene rings directly connected.
    biphenyl_smarts = "c1ccc(cc1)-c2ccccc2"
    biphenyl_core = Chem.MolFromSmarts(biphenyl_smarts)
    
    cores = [
        ("dibenzodioxin", dioxin_core),
        ("dibenzofuran", dibenzofuran_core),
        ("biphenyl", biphenyl_core)
    ]
    
    # Try to find one core that meets our halogenation threshold.
    for core_name, core_pattern in cores:
        if core_pattern is None:
            continue
        matches = mol.GetSubstructMatches(core_pattern)
        if matches:
            # For each match (which is a tuple of atom indices forming the core)
            for match in matches:
                # Create a set for quick lookup.
                core_atom_idxs = set(match)
                # Count halogen atoms attached to any atom in the core.
                attached_halogen_count = 0
                for idx in match:
                    atom = mol.GetAtomWithIdx(idx)
                    for nbr in atom.GetNeighbors():
                        # Count the neighbor if it is a halogen and not part of the core.
                        if nbr.GetIdx() not in core_atom_idxs and nbr.GetAtomicNum() in halogen_znums:
                            attached_halogen_count += 1
                if attached_halogen_count >= 2:
                    return True, (f"CORRECT: Contains a {core_name} core with "
                                  f"{attached_halogen_count} direct halogen substituents")
    # If no match qualifies then provide a reason.
    # Also, if no core was found at all.
    core_found = any(mol.HasSubstructMatch(core) for (_, core) in cores if core is not None)
    if not core_found:
        return False, "Does not contain a dioxin, dibenzofuran, or biphenyl core"
    return False, "Core found but not enough halogen substituents directly attached to the aromatic core"

# (Optional) For testing purposes:
if __name__ == "__main__":
    test_examples = [
        ("1,2,3,7,8-Pentachlorodibenzodioxin", "Clc1cc2Oc3cc(Cl)c(Cl)c(Cl)c3Oc2cc1Cl"),
        ("2,2',3,3',5,5'-hexachlorobiphenyl", "Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)cc(Cl)c1Cl"),
        ("Formicamycin I", "ClC1=C(O)C(Cl)=C2C([C@H]3CC=4C(Cl)=C(OC)C=C(C4C([C@]3(C(C2=C1O)=O)O)=O)C5=C(O)C=C(OC)C(=C5C)Cl)(C)C"),
        # Below are some that were false negatives previously.
        ("3,5-dichlorobiphenyl", "Clc1cc(Cl)cc(c1)-c1ccccc1"),
        ("niclofolan", "OC1=C(C=C(C=C1[N+]([O-])=O)Cl)C2=C(O)C([N+](=O)[O-])=CC(=C2)Cl"),
        ("4,4'-dichlorobiphenyl", "Clc1ccc(cc1)-c1ccc(Cl)cc1"),
        # A false positive example from before (a larger molecule) that should now be rejected.
        ("Chlorobiphenyl-desleucyl-vancomycin", "C[C@@H]1O[C@H](C[C@](C)(NCc2ccc(cc2)-c2ccc(Cl)cc2)[C@@H]1O)O[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1Oc1c2Oc3ccc(cc3Cl)[C@@H](O)[C@@H](N)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H]3c(c2)cc1Oc1ccc(cc1Cl)[C@@H](O)[C@@H]1NC(=O)[C@H](NC3=O)c2ccc(O)c(c2)-c2c(O)cc(O)cc2[C@H](NC1=O)C(O)=O")
    ]
    for name, smi in test_examples:
        result, reason = is_polychlorinated_dibenzodioxines_and_related_compounds(smi)
        print(f"{name}: {result} ({reason})")