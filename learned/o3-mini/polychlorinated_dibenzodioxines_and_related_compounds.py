"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Polychlorinated dibenzodioxines and related compounds.
These include organochlorine/organobromine compounds such as 
polychlorinated dibenzodioxins, dibenzofurans, and polyhalogenated biphenyls.
Improvements over a previous attempt:
  • For each aromatic core match, we now require that the number of atoms in the match is exactly the number of atoms
    defined by the core SMARTS – this helps ensure that the core is not simply a fragment of a larger aromatic network.
  • We also check that none of the core atoms has an extra directly attached aromatic substituent
    (unless that substituent is a halogen).
  • Each core type requires a minimal number of attached halogen substituents:
        – dibenzodioxin:  threshold >= 3
        – dibenzofuran: threshold >= 3
        – biphenyl:     threshold >= 2
"""

from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to polychlorinated dibenzodioxines and related compounds.
    
    The procedure is:
      1. Parse and sanitize the molecule.
      2. For each candidate aromatic core substructure (dibenzodioxin, dibenzofuran, biphenyl):
         a. Get all substructure matches.
         b. For each match, require that the number of atoms in the match exactly equals the number of atoms 
            specified by the SMARTS pattern (to avoid matching only part of an extended aromatic system).
         c. Check that none of the core atoms has an additional aromatic substituent that is not itself a halogen.
         d. Count the number of directly attached halogen substituents.
         e. If the number of attached halogens is at or above the core‐specific threshold, accept the match.
      3. If no match passes, return a negative classification.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule belongs to this chemical class; False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize molecule to ensure aromaticity and valence assignment.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"
    
    # Define the set of halogens (atomic numbers) that are acceptable substituents.
    halogen_znums = {17, 35}  # Cl and Br

    # Define candidate aromatic cores: each tuple is (core_name, SMARTS pattern, minimum halogen threshold)
    cores = [
        ("dibenzodioxin", "c1ccc2c(c1)Oc3ccccc3O2", 3),
        ("dibenzofuran", "c1ccc2occc2c1", 3),
        ("biphenyl",    "c1ccc(cc1)-c2ccccc2", 2)
    ]
    
    # For each core, look for matches in the molecule
    for core_name, smarts, threshold in cores:
        core_pattern = Chem.MolFromSmarts(smarts)
        if core_pattern is None:
            continue
        # Get all substructure matches
        matches = mol.GetSubstructMatches(core_pattern)
        for match in matches:
            # Only consider matches that cover the entire core (avoid partial matches from extended networks)
            if len(match) != core_pattern.GetNumAtoms():
                continue
            
            core_atom_idxs = set(match)
            # Check for any extra aromatic substituent attached to the core atoms.
            # We disallow any neighbor (outside the core) which is aromatic and is not a halogen.
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
                # This match is not a "clean core", so try the next match.
                continue
            
            # Count the number of halogen substituents directly attached to the atoms in the core.
            attached_halogen_count = 0
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in core_atom_idxs and nbr.GetAtomicNum() in halogen_znums:
                        attached_halogen_count += 1
            
            # If enough halogen substituents are found for this core, classify as a match.
            if attached_halogen_count >= threshold:
                return True, (f"CORRECT: Contains a {core_name} core with "
                              f"{attached_halogen_count} direct halogen substituents")
    
    # If we have found no match that passes the criteria, signal that either no core was found
    # or the cores were extended (had extra aromatic substituents) or did not have enough halogens.
    core_found = any(mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)) for (_, smarts, _) in cores)
    if not core_found:
        return False, "Does not contain a dibenzodioxin, dibenzofuran, or biphenyl core"
    
    return False, ("Core found but either the core is extended by additional aromatic substituents "
                   "or does not have enough direct halogen substituents")
    
# (Optional) Testing examples
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
        # An example of a compound that should not be classified:
        ("GSK805", "O=C(NC1=CC(Cl)=C(C2=CC=CC=C2OC(F)(F)F)C(Cl)=C1)CC3=CC=C(S(=O)(CC)=O)C=C3")
    ]
    
    for name, smi in test_examples:
        result, reason = is_polychlorinated_dibenzodioxines_and_related_compounds(smi)
        print(f"{name}: {result} ({reason})")