"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: Mineral nutrient
Definition: A mineral nutrient is defined as an inorganic nutrient that must be ingested 
and absorbed in adequate amounts to satisfy a wide range of essential metabolic and/or 
structural functions in the human body. In our heuristic, in addition to requiring one of 
a set of known metal (or metalloid) ions, we reject molecules that appear “organic” or 
coordination-complex by checking if the compound contains rings or triple bonds. We also 
add a special rule to avoid classifying simple Al–chloride compounds as nutrient.
"""

from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines whether a molecule (given as SMILES) is a mineral nutrient.
    Our heuristic is as follows:
      • The molecule must be parsed correctly and contain at least one metal from a predefined set.
      • We expect a mineral nutrient to be an inorganic salt or simple ion pairing.
      • Therefore, if the molecule contains any rings (cyclic substructures) or triple bonds 
        (as found in some coordination complexes) we reject it.
      • We also reject molecules containing iodine (I) since none of our nutrient examples have I.
      • Finally, we add a specific rule: although aluminium sulfate (and similar) is allowed,
        a simple aluminium–chloride such as AlCl3 is rejected.
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) where the bool indicates whether the molecule is classified as a mineral nutrient,
               and the string provides the reason.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a set of allowed metal elements (based on our examples)
    allowed_metals = {"Fe", "Zn", "Ca", "K", "Mg", "Na", "Ba", "Al", "Cs", "Pd", "Sb", "La"}
    found_metals = set()
    nonmetal_set = set()
    
    # For every atom in the molecule, record its element and whether it is one of our allowed metals.
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in allowed_metals:
            found_metals.add(symbol)
        else:
            nonmetal_set.add(symbol)
    
    # If no allowed metal is present, this is not a mineral nutrient.
    if not found_metals:
        return False, "No recognized mineral metal found; not a mineral nutrient"
    
    # Check if the molecule contains any rings.
    # (Many complex organic/coordination compounds contain rings. In our examples, the nutrient salts do not.)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s); appears to be a complex organic or coordination compound"
    
    # Check for any triple bonds. (Many false positives featured C#N triple bonds.)
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            return False, "Molecule contains triple bonds; appears to be a complex coordination compound"
    
    # If any iodine atoms are present (besides metals) we reject the compound.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "I":
            return False, "Contains iodine; not typical for a mineral nutrient"
    
    # Special rule for aluminium:
    # There are examples where aluminium in a sulfate salt is nutrient, but AlCl3 type compounds are not.
    # Here, if Al is the only metal and the only non-metal elements are Cl, we reject.
    if "Al" in found_metals:
        # Using set comparison: if the nonmetal set is exactly {"Cl"} then this is probably just AlCl3 (or similar)
        if nonmetal_set == {"Cl"}:
            return False, "Aluminium present with only chloride as counterion; not considered a mineral nutrient"
    
    # Also, if there are no non-metal atoms at all (for example, a pure metal cluster) reject.
    heavy_nonmetal_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() not in allowed_metals]
    if len(heavy_nonmetal_atoms) == 0:
        return False, "No counter ion (non-metal atoms) found; not a typical salt for a mineral nutrient"
    
    # If we made it here the molecule appears to be a simple salt (or ion‐pair) with allowed metal(s) and
    # without features marking it as an organic or coordination complex.
    return True, "Contains mineral metal(s): " + ", ".join(sorted(found_metals))

# Example test calls:
if __name__ == "__main__":
    test_smiles = [
        "[Fe+3].[O-]P([O-])(=O)[O-]",  # iron(3+) phosphate (should be True)
        "Cl[Al](Cl)Cl",               # aluminium trichloride (should be False)
        "[Na+].[Na+].O=N[Fe--](C#N)(C#N)(C#N)(C#N)C#N",  # sodium nitroprusside (should be False)
        "[Mg++].CCCCCCCCCCCCCCCCCC([O-])=O.CCCCCCCCCCCCCCCCCC([O-])=O",  # magnesium distearate (True)
    ]
    for s in test_smiles:
        res, reason = is_mineral_nutrient(s)
        print(f"SMILES: {s}\nResult: {res} \nReason: {reason}\n{'-'*40}")