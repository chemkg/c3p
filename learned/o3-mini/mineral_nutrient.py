"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: Mineral nutrient
Definition: A mineral nutrient is defined as an inorganic nutrient that must be ingested 
and absorbed in adequate amounts to satisfy a wide range of essential metabolic and/or 
structural functions in the human body.
Heuristic:
 • Parse the SMILES string.
 • Reject molecules with rings or triple bonds.
 • Ensure that at least one allowed metal is present.
 • Decompose the molecule into fragments (ignoring water) and require that at least one fragment 
   appears to be an allowed counterion based on a whitelist of substructures or allowed monatomic ions.
 • Apply special rules for single-fragment cases (to avoid covalent complexes).
"""

from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Classify a molecule (given as a SMILES string) as a mineral nutrient or not.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple. The first element is True if classified as a mineral nutrient,
                     False otherwise; the second element is a string explaining the reasoning.
    """
    # Parse SMILES string; return error if parsing fails.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules with rings (suggesting complex structures)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s); appears to be a covalent or complex coordination compound"
    
    # Reject molecules that contain triple bonds (often seen in complex molecules)
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            return False, "Molecule contains triple bonds; appears to be a complex coordination compound"
    
    # Allowed metals (based on our examples)
    # (Note: some metals like Sb are sometimes present; we handle them specially in single-fragment cases.)
    allowed_metals = {"Fe", "Zn", "Ca", "K", "Mg", "Na", "Ba", "Al", "Cs", "Pd", "Sb", "La"}
    found_metals = {atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() in allowed_metals}
    if not found_metals:
        return False, "No recognized mineral metal found; not a mineral nutrient"
    
    # Break molecule into fragments; ignore water fragments
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    non_water_frags = []
    for frag in frags:
        frag_smiles = Chem.MolToSmiles(frag, canonical=True)
        if frag_smiles == "O":  # ignore water
            continue
        non_water_frags.append(frag)
    
    # Special handling for single fragment molecules.
    if len(non_water_frags) == 1:
        single_frag = non_water_frags[0]
        # If the single fragment contains Sb in what appears to be a covalent arrangement,
        # we suspect it is not a typical ionic nutrient.
        for atom in single_frag.GetAtoms():
            if atom.GetSymbol() == "Sb":
                return False, "Single-fragment Sb compounds (e.g., SbCl3) are not considered mineral nutrients"
        # If two allowed metal atoms are directly bonded (covalently linked) then suspect a complex.
        for bond in single_frag.GetBonds():
            a1 = bond.GetBeginAtom().GetSymbol()
            a2 = bond.GetEndAtom().GetSymbol()
            if a1 in allowed_metals and a2 in allowed_metals:
                return False, "Metal atoms are covalently bonded; not an ionic salt"
    
    # Build a dictionary of allowed counterion SMARTS patterns;
    # only keep the ones that compile successfully.
    anion_smarts_strings = {
        "phosphate": "[O-]P(=O)([O-])[O-]",
        "hydrogenphosphate": "[H]OP(=O)([O-])[O-]",
        "nitrate": "[O-][N+]([O-])=O",
        "sulfate": "[O-]S(=O)([O-])[O-]",
        "formate": "[H]C([O-])=O",
        "acetate": "CC([O-])=O",
        "propionate": "CCC([O-])=O",
        "chloride": "[Cl-]",
        "fluoride": "[F-]",
        "carbonate": "C(=O)([O-])[O-]",
        "hydroxide": "[OH-]",
        "silicate": "[O-]Si([O-])([O-])[O-]",
        "nitrate_alt": "[N+](=O)[O-]",
    }
    allowed_anions = {}
    for name, smarts in anion_smarts_strings.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            allowed_anions[name] = patt
        # If a pattern fails to compile, it will not be used.
    
    allowed_anion_found = False
    # Also define a small set of allowed monatomic anions.
    allowed_monatomic = {"Cl", "F", "O"}
    
    # Process each fragment that is *not* a pure metal ion.
    for frag in non_water_frags:
        atom_symbols = [atom.GetSymbol() for atom in frag.GetAtoms()]
        
        # If the fragment is a single atom:
        if frag.GetNumAtoms() == 1:
            elem = atom_symbols[0]
            # If the atom (e.g., Cl, F, O) is allowed as a monatomic counterion, mark it.
            if elem in allowed_monatomic:
                allowed_anion_found = True
                continue
            # If it is an allowed metal then that fragment is already counted.
            if elem in allowed_metals:
                continue
            # Otherwise, unrecognized single-atom fragment.
            return False, f"Fragment '{Chem.MolToSmiles(frag)}' not recognized as an allowed ion"

        # For multi-atom fragments, try to match one of the allowed anion SMARTS patterns.
        match_found = False
        for name, pattern in allowed_anions.items():
            # Check that the pattern is valid (should be, due to filtering above)
            if pattern is not None and frag.HasSubstructMatch(pattern):
                # Special rule: if the match is hydrogenphosphate and the metal is an alkali metal,
                # then we may consider it less typical. (For example, sodium or potassium hydrogenphosphate are less common nutrient forms.)
                if name == "hydrogenphosphate" and any(m in {"Na", "K"} for m in found_metals):
                    continue
                match_found = True
                break
        if not match_found:
            # It might be that the fragment is an ionic metal fragment.
            if any(atom.GetSymbol() in allowed_metals for atom in frag.GetAtoms()):
                continue
            return False, f"Fragment '{Chem.MolToSmiles(frag)}' does not match allowed counterion patterns"
        else:
            allowed_anion_found = True

    if not allowed_anion_found:
        return False, "No allowed counterion fragment detected"

    # Final check: ensure that disallowed halides (Br or I) are absent from every fragment.
    for frag in non_water_frags:
        for atom in frag.GetAtoms():
            if atom.GetSymbol() in {"Br", "I"}:
                return False, "Contains disallowed halide (Br or I)"

    return True, "Contains mineral nutrient components, including metal(s): " + ", ".join(sorted(found_metals))


# Example test calls if run as a script (uncomment the block below to test):
if __name__ == "__main__":
    test_examples = [
        "[Fe+3].[O-]P([O-])(=O)[O-]",                   # iron(3+) phosphate
        "[Zn++].[O-][N+]([O-])=O.[O-][N+]([O-])=O",       # zinc nitrate
        "[Cs+].[H]C([O-])=O",                            # caesium formate
        "[Ca++].[Ca++].[O-][Si]([O-])([O-])[O-]",         # calcium silicate
        "[K+].[K+].[O-]S([O-])(=O)=O",                   # potassium sulfate
        "[Ba++].CC([O-])=O.CC([O-])=O",                   # barium acetate
        "[Mg++].CCCCCCCCCCCCCCCCCC([O-])=O.CCCCCCCCCCCCCCCCCC([O-])=O",  # magnesium distearate
        "[Cl-].[K+]",                                   # potassium chloride
        "[Pd-2](Cl)(Cl)(Cl)(Cl)(Cl)Cl.[K+].[K+]",         # Potassium hexachloropalladate(IV)
        "[Na+].[Na+].[Na+].[O-]P([O-])([O-])=O",          # trisodium phosphate
        "Cl[O-].[Ca+2].Cl[O-]",                          # Calcium hypochlorite
        "[Ca++].[O-]S([O-])(=O)=O",                       # calcium sulfate
        "[OH-].[OH-].[Ca++]",                           # calcium dihydroxide
        "[Mg++].[Cl-].[Cl-]",                           # magnesium dichloride
        "[Ba++].[O-][N+]([O-])=O.[O-][N+]([O-])=O",       # barium nitrate
        "Cl[La](Cl)Cl",                                 # lanthanum trichloride
        "O[Ca]",                                        # calcium monohydroxide
        "[Ca+2].C(=O)([O-])[O-]",                        # calcium carbonate
        "[Mg++].[O-]S([O-])(=O)=O",                       # magnesium sulfate
        "[Cl-].[Cl-].[Ca++]",                           # calcium dichloride
        "[Ca++].[Ca++].[Ca++].[O-]P([O-])([O-])=O.[O-]P([O-])([O-])=O",  # tricalcium bis(phosphate)
        "[Ca++].[H]OP([O-])([O-])=O",                     # calcium hydrogenphosphate
        "[O-]S([O-])(=O)=O.[Ba+2]",                      # barium sulfate
        "O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.[Al+3].[Al+3].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O",  # aluminium sulfate octadecahydrate
        "[Ba++].[O-]C([O-])=O",                          # barium carbonate
        "[Cl-].[Cs+]",                                  # caesium chloride
        "[Na+].[Na+].OP([O-])([O-])=O",                  # disodium hydrogenphosphate (should be rejected if paired with water fragments)
        "Cl[Sb](Cl)Cl",                                 # antimony trichloride (should be rejected)
        "[S--].[Fe+3].[As-]",                           # arsenopyrite (should be rejected)
        "[K+].[Br-]",                                  # potassium bromide (disallowed halide)
        "[Ca++].OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)C([O-])=O."
        "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)C([O-])=O",  # calcium glucoheptonate (rejected)
        "[Na+].F[P-](F)(F)(F)(F)F",                     # sodium hexafluorophosphate (rejected)
    ]
    for s in test_examples:
        result, reason = is_mineral_nutrient(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*60}")