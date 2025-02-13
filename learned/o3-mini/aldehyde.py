"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: Aldehyde (RC(=O)H)
Definition: A compound in which a carbonyl group (C=O) is bonded to exactly one hydrogen atom and one R group.
This function will not consider formaldehyde (H2C=O) or carboxylic acids (where the noncarbonyl substituent is oxygen) as aldehydes.
"""

from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    For this purpose, an aldehyde is defined as a molecule having a formyl group – a carbonyl (C=O) 
    where the carbon has exactly one hydrogen and exactly one heavy-atom substituent (R group),
    with the heavy substituent being nonoxygen (to avoid acids, esters, etc). 
     
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an aldehyde, False otherwise.
        str: Explanation for the classification.
    """

    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Add explicit hydrogens so that hydrogen counts are correct.
    mol = Chem.AddHs(mol)

    # Loop over all atoms and look for a carbon atom that could be the aldehyde carbon.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue

        # Find the double-bonded oxygen neighbors (the carbonyl oxygen).
        carbonyl_oxygens = []
        for bond in atom.GetBonds():
            # Check for a double bond.
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetSymbol() == "O":
                    carbonyl_oxygens.append(nbr)
        # For a valid aldehyde, we require exactly one carbonyl oxygen.
        if len(carbonyl_oxygens) != 1:
            continue

        # Now consider the other neighbors (excluding the carbonyl oxygen).
        # They should consist of exactly one heavy atom (the R group) and the hydrogens.
        heavy_neighbors = []
        for nbr in atom.GetNeighbors():
            # Skip the carbonyl oxygen; we already accounted for it.
            if nbr.GetIdx() in [o.GetIdx() for o in carbonyl_oxygens]:
                continue
            # Only count heavy atoms (atomic number > 1).
            if nbr.GetAtomicNum() > 1:
                heavy_neighbors.append(nbr)

        # According to the definition, the carbon in RC(=O)H should have exactly one heavy neighbor.
        if len(heavy_neighbors) != 1:
            continue

        # Exclude cases where the noncarbonyl substituent is oxygen (e.g. carboxylic acid groups).
        if heavy_neighbors[0].GetSymbol() == "O":
            continue

        # Lastly, check that the carbon has exactly one hydrogen attached.
        # GetTotalNumHs counts both implicit and explicit hydrogens.
        if atom.GetTotalNumHs() != 1:
            continue

        # If all conditions are met, we have detected an aldehyde group.
        return True, ("Aldehyde group detected: the carbonyl carbon (C=O) is bonded to exactly one hydrogen "
                      "and one non-oxygen substituent (R group).")

    # No atom meeting the aldehyde criteria was found.
    return False, "No aldehyde functional group (RC(=O)H) detected."

# Example usage:
if __name__ == '__main__':
    # Test a few sample SMILES strings from the provided list.
    smiles_examples = [
        "O=CC1=C(O)C=C(C#CC(=C)C)C(=C1)O",         # 2,5-dihydroxy-4-(3-methylbut-3-en-1-ynyl)benzaldehyde
        "C[C@H](C=O)C(O)=O",                         # (2R)-2-methyl-3-oxopropanoic acid
        "CN(CCC=O)N=O",                              # 3-(N-Nitrosomethylamino)propionaldehyde
        "O=CCCCCCCCCC/C=C/C=C\\CC",                   # 11E,13Z-Hexadecadienal
        "[H]C(=O)CC#C",                              # but-3-ynal
        "[H]C(=O)C(C)CCCCCCCCCCCCC",                  # 2-methylpentadecanal
        "[H]C(=O)c1ccc2ccccc2c1",                     # 2-naphthaldehyde
        "[H]C(=O)C(C)C",                             # A simple aldehyde example (isobutyraldehyde would be RC(=O)H)
    ]
    for s in smiles_examples:
        result, reason = is_aldehyde(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*60}")