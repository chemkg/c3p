"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Polychlorinated dibenzodioxines and related compounds.
These include organochlorine/organobromine compounds such as 
polychlorinated dibenzodioxins, dibenzofurans, and polyhalogenated biphenyls.
Improvements over a previous attempt:
  - For each aromatic core match, we now check that none of the core atoms have an extra aromatic substituent 
    (other than a halogen). This “clean core” criterion helps rule out large molecules in which a biphenyl-like substructure
    is part of an extended aromatic network.
  - We set core‐dependent thresholds for the number of halogen substituents:
      • dibenzodioxin: >= 3 attached halogens
      • dibenzofuran:  >= 3 attached halogens
      • biphenyl:     >= 2 attached halogens
"""

from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to polychlorinated dibenzodioxines and related compounds.
    
    The classification procedure is:
      1. Parse the molecule and sanitize it.
      2. Identify one of three aromatic cores (dibenzodioxin/dibenzofuran/biphenyl) by substructure matching.
      3. For each substructure match (core):
           a. Check that none of the core atoms is substituted by an extra aromatic fragment (i.e. a neighbor that is aromatic but not a halogen).
           b. Count the number of halogen substituents (Cl or Br) attached directly to the core.
      4. If a “clean” core is found with enough attached halogens (thresholds per core type), the function returns True.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule belongs to this chemical class; False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure aromaticity and proper valence assignment.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"
    
    # Set of halogens of interest (atomic numbers)
    halogen_znums = {17, 35}  # Cl and Br

    # Define the target aromatic cores with SMARTS and the minimal halogen thresholds.
    cores = [
        # (core_name, SMARTS pattern, minimum required attached halogens)
        ("dibenzodioxin", "c1ccc2c(c1)Oc3ccccc3O2", 3),
        ("dibenzofuran", "c1ccc2occc2c1", 3),
        ("biphenyl",    "c1ccc(cc1)-c2ccccc2", 2)
    ]
    
    # For each core type, attempt to match
    for core_name, smarts, threshold in cores:
        core_pattern = Chem.MolFromSmarts(smarts)
        if core_pattern is None:
            # Skip if pattern generation failed.
            continue
        matches = mol.GetSubstructMatches(core_pattern)
        if matches:
            # Evaluate each match independently.
            for match in matches:
                core_atom_idxs = set(match)
                # Check for any "extra" aromatic substitution attached directly to a core atom.
                # We disallow neighbors (not part of the match) that are aromatic and NOT a halogen.
                clean_core = True
                for idx in match:
                    atom = mol.GetAtomWithIdx(idx)
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() not in core_atom_idxs:
                            if nbr.GetIsAromatic() and nbr.GetAtomicNum() not in halogen_znums:
                                clean_core = False
                                break
                    if not clean_core:
                        break
                if not clean_core:
                    # Skip this match if the core is extended by extra aromatic fragments.
                    continue
                
                # Count the number of halogen substituents directly attached to the core.
                attached_halogen_count = 0
                for idx in match:
                    atom = mol.GetAtomWithIdx(idx)
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() not in core_atom_idxs and nbr.GetAtomicNum() in halogen_znums:
                            attached_halogen_count += 1
                if attached_halogen_count >= threshold:
                    return True, (f"CORRECT: Contains a {core_name} core with "
                                  f"{attached_halogen_count} direct halogen substituents")
    
    # If no qualifying core is found, check if any core was found at all:
    core_found = any(mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)) for (_, smarts, _) in cores)
    if not core_found:
        return False, "Does not contain a dibenzodioxin, dibenzofuran, or biphenyl core"
    
    return False, ("Core found but either the core is extended by additional aromatic substituents "
                   "or does not have enough direct halogen substituents")
    
# (Optional) For testing purposes.
if __name__ == "__main__":
    test_examples = [
        ("1,2,3,7,8-Pentachlorodibenzodioxin", "Clc1cc2Oc3cc(Cl)c(Cl)c(Cl)c3Oc2cc1Cl"),
        ("2,2',3,3',5,5'-hexachlorobiphenyl", "Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)cc(Cl)c1Cl"),
        ("Formicamycin I", "ClC1=C(O)C(Cl)=C2C([C@H]3CC=4C(Cl)=C(OC)C=C(C4C([C@]3(C(C2=C1O)=O)O)=O)C5=C(O)C=C(OC)C(=C5C)Cl)(C)C"),
        ("PCB138", "ClC1=C(C=2C(Cl)=CC(Cl)=C(Cl)C2)C=CC(Cl)=C1Cl"),
        ("1,2,3,4,7,8,9-Heptachlorodibenzofuran", "Clc1cc2oc3c(Cl)c(Cl)c(Cl)c(Cl)c3c2c(Cl)c1Cl"),
        ("2,3',4,4',5-Pentachlorobiphenyl", "Clc1ccc(cc1Cl)-c1cc(Cl)c(Cl)cc1Cl"),
        ("2,2',4,4',5,5'-hexabromobiphenyl", "C1(C2=C(C=C(C(=C2)Br)Br)Br)=C(C=C(C(=C1)Br)Br)Br"),
        ("3,5-dichlorobiphenyl", "Clc1cc(Cl)cc(c1)-c1ccccc1"),
        ("niclofolan", "OC1=C(C=C(C=C1[N+]([O-])=O)Cl)C2=C(O)C([N+](=O)[O-])=CC(=C2)Cl"),
        ("4,4'-dichlorobiphenyl", "Clc1ccc(cc1)-c1ccc(Cl)cc1"),
        # False positive candidate from previous attempt:
        ("GSK805", "O=C(NC1=CC(Cl)=C(C2=CC=CC=C2OC(F)(F)F)C(Cl)=C1)CC3=CC=C(S(=O)(CC)=O)C=C3")
    ]
    
    for name, smi in test_examples:
        result, reason = is_polychlorinated_dibenzodioxines_and_related_compounds(smi)
        print(f"{name}: {result} ({reason})")