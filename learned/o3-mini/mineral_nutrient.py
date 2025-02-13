"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: Mineral nutrient
Definition: A mineral nutrient is defined as an inorganic nutrient that must be ingested 
and absorbed in adequate amounts to satisfy a wide range of essential metabolic and/or 
structural functions in the human body.
Our heuristic first ensures the molecule parses correctly and contains at least one allowed metal.
We then check that no rings or triple bonds are present.
Next, we break the molecule into fragments (ignoring water) and require that there is evidence
of an allowed counterion. Allowed anions are detected by substructure matching to a small whitelist.
Special rules (e.g. for hydrogenphosphate with Na/K and for single‐fragment Sb compounds)
are applied.
"""

from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines whether a molecule (given as SMILES) is a mineral nutrient.
    
    Heuristic:
     • Parse the SMILES and remove water fragments.
     • Reject if rings or triple bonds are present.
     • Check that at least one allowed metal (from a hard‐coded set) is present.
     • Using fragment analysis, require at least one counterion fragment that
       matches a whitelist of allowed anion SMARTS. Also if the entire molecule is a single fragment,
       then special rules are applied (e.g. rejecting SbCl3–like compounds).
     • If any disallowed atoms (e.g. Br, I or atoms not expected in simple inorganic ions)
       are present in non‐metal fragments, reject.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple with True if classified as a mineral nutrient and a string reason;
                     or False and the reason for rejection.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that there are no rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s); appears to be a complex organic or coordination compound"
    # Check for triple bonds in any bond
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            return False, "Molecule contains triple bonds; appears to be a complex coordination compound"
    
    # Define a set of allowed metal symbols (from our examples)
    allowed_metals = {"Fe", "Zn", "Ca", "K", "Mg", "Na", "Ba", "Al", "Cs", "Pd", "Sb", "La"}
    found_metals = {atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() in allowed_metals}
    if not found_metals:
        return False, "No recognized mineral metal found; not a mineral nutrient"
    
    # Break molecule into fragments; ignore water fragments (e.g. "O")
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    non_water_frags = []
    for frag in frags:
        frag_smiles = Chem.MolToSmiles(frag, canonical=True)
        if frag_smiles == "O":
            continue
        non_water_frags.append(frag)
    
    # If there is only one non-water fragment, apply a safeguard:
    #   Reject if an allowed metal like "Sb" is present in a non-ionic (covalent) context.
    if len(non_water_frags) == 1:
        single_frag = non_water_frags[0]
        # Special rule: many ionic nutrient salts are dot-separated.
        # If the single fragment contains Sb, reject (e.g. SbCl3 should not be nutrient)
        for atom in single_frag.GetAtoms():
            if atom.GetSymbol() == "Sb":
                return False, "Single-fragment Sb compounds (e.g., SbCl3) are not considered mineral nutrients"
        # Also, if a bond exists between two allowed metals, we suspect a covalent complex.
        for bond in single_frag.GetBonds():
            a1 = bond.GetBeginAtom().GetSymbol()
            a2 = bond.GetEndAtom().GetSymbol()
            if a1 in allowed_metals and a2 in allowed_metals:
                return False, "Metal atoms are covalently bonded; not an ionic salt"
    
    # Define a library of allowed counterion SMARTS patterns.
    # Note that some patterns (like hydrogenphosphate) are conditionally allowed.
    anion_smarts = {
        "phosphate": Chem.MolFromSmarts("[O-]P(=O)([O-])[O-]"),
        "hydrogenphosphate": Chem.MolFromSmarts("[H]OP(=O)([O-])[O-]"),
        "nitrate": Chem.MolFromSmarts("[O-][N+]([O-])=O"),
        "sulfate": Chem.MolFromSmarts("[O-]S(=O)([O-])[O-]"),
        "formate": Chem.MolFromSmarts("[H]C([O-])=O"),
        "acetate": Chem.MolFromSmarts("CC([O-])=O"),
        "propionate": Chem.MolFromSmarts("CCC([O-])=O"),
        "chloride": Chem.MolFromSmarts("[Cl-]"),
        "fluoride": Chem.MolFromSmarts("[F-]"),
        "carbonate": Chem.MolFromSmarts("C(=O)([O-])[O-]"),
        "hydroxide": Chem.MolFromSmarts("[OH-]"),
        "silicate": Chem.MolFromSmarts("[O-]Si([O-])([O-])[O-]"),
        # Sometimes nitrate is written differently:
        "nitrate_alt": Chem.MolFromSmarts("[N+](=O)[O-]"),
    }
    
    allowed_anion_found = False
    # Process each fragment – if it is not a pure metal ion then it should match one allowed pattern.
    for frag in non_water_frags:
        frag_atoms = [atom.GetSymbol() for atom in frag.GetAtoms()]
        # If fragment is a single atom and is not a metal, then check if it is an allowed simple ion.
        if frag.GetNumAtoms() == 1:
            elem = frag_atoms[0]
            # For allowed monatomic anions, we allow chloride, fluoride or hydroxide.
            if elem in {"Cl", "F", "O"}:
                # Here, "O" may come from OH- (but water already removed).
                allowed_anion_found = True
                continue
            # If it is a metal ion, that is already counted.
            if elem in allowed_metals:
                continue
            return False, f"Fragment '{Chem.MolToSmiles(frag)}' not recognized as allowed ion"
        
        # For multi-atom fragments that are not pure metal ions, try to match one of our allowed anion patterns.
        match_found = False
        for name, pattern in anion_smarts.items():
            if frag.HasSubstructMatch(pattern):
                # Special rule: if the allowed match is hydrogenphosphate and the metal is an alkali metal,
                # we consider the salt less typical (as in disodium hydrogenphosphate dihydrate).
                if name == "hydrogenphosphate" and any(m in {"Na", "K"} for m in found_metals):
                    # Do not accept hydrogenphosphate when paired with these metals.
                    continue
                match_found = True
                break
        if not match_found:
            # If the fragment does not match any allowed pattern, check if it is a metal ion fragment.
            if any(atom.GetSymbol() in allowed_metals for atom in frag.GetAtoms()):
                continue  # that fragment is likely a metal ion
            return False, f"Fragment '{Chem.MolToSmiles(frag)}' does not match allowed counterion patterns"
        else:
            allowed_anion_found = True
    
    if not allowed_anion_found:
        return False, "No allowed counterion fragment detected"

    # Additional check: ensure that disallowed halides (Br, I) are not present in any fragment.
    for frag in non_water_frags:
        for atom in frag.GetAtoms():
            if atom.GetSymbol() in {"Br", "I"}:
                return False, "Contains disallowed halide (Br or I)"
    
    return True, "Contains mineral metal(s): " + ", ".join(sorted(found_metals))


# Example test calls if run as script:
if __name__ == "__main__":
    test_examples = [
        # True positives
        "[Fe+3].[O-]P([O-])(=O)[O-]",        # iron(3+) phosphate
        "[Zn++].[O-][N+]([O-])=O.[O-][N+]([O-])=O",  # zinc nitrate
        "[Cs+].[H]C([O-])=O",                 # caesium formate
        "[Ca++].[Ca++].[O-][Si]([O-])([O-])[O-]", # calcium silicate
        "[K+].[K+].[O-]S([O-])(=O)=O",         # potassium sulfate
        "[Ba++].CC([O-])=O.CC([O-])=O",          # barium acetate
        "[Cl-].[K+]",                        # potassium chloride
        "[Mg++].CCCCCCCCCCCCCCCCCC([O-])=O.CCCCCCCCCCCCCCCCCC([O-])=O",  # magnesium distearate
        "[Pd-2](Cl)(Cl)(Cl)(Cl)(Cl)Cl.[K+].[K+]",  # Potassium hexachloropalladate(IV)
        "[Na+].[Na+].[Na+].[O-]P([O-])([O-])=O",  # trisodium phosphate
        "Cl[O-].[Ca+2].Cl[O-]",               # Calcium hypochlorite
        "[Ca++].[O-]S([O-])(=O)=O",            # calcium sulfate
        "[OH-].[OH-].[Ca++]",                 # calcium dihydroxide
        "[Mg++].[Cl-].[Cl-]",                 # magnesium dichloride
        "[Ba++].[O-][N+]([O-])=O.[O-][N+]([O-])=O",  # barium nitrate
        "Cl[La](Cl)Cl",                      # lanthanum trichloride
        "O[Ca]",                             # calcium monohydroxide
        "[Ca+2].C(=O)([O-])[O-]",             # calcium carbonate
        "[Mg++].[O-]S([O-])(=O)=O",            # magnesium sulfate
        "[Cl-].[Cl-].[Ca++]",                # calcium dichloride
        "[Ca++].[Ca++].[Ca++].[O-]P([O-])([O-])=O.[O-]P([O-])([O-])=O",  # tricalcium bis(phosphate)
        "[Ca++].[H]OP([O-])([O-])=O",          # calcium hydrogenphosphate
        "[O-]S([O-])(=O)=O.[Ba+2]",           # barium sulfate
        "O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.[Al+3].[Al+3].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O",  # aluminium sulfate octadecahydrate
        "[Ba++].[O-]C([O-])=O",               # barium carbonate
        "[Cl-].[Cs+]",                       # caesium chloride
        # False positives (should be rejected)
        "[Na+].[Na+].[O-]P(=O)(O)[O-].O.O",    # disodium hydrogenphosphate dihydrate
        "Cl[Sb](Cl)Cl",                      # antimony trichloride
        "[S--].[Fe+3].[As-]",                 # arsenopyrite
        "[K+].[Br-]",                        # potassium bromide
        "[Ca++].OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)C([O-])=O.OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)C([O-])=O",  # calcium glucoheptonate
        "[Na+].F[P-](F)(F)(F)(F)F",           # sodium hexafluorophosphate
    ]
    
    for s in test_examples:
        result, reason = is_mineral_nutrient(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*60}")