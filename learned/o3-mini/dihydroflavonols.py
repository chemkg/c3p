"""
Classifies: CHEBI:48039 dihydroflavonols
"""
#!/usr/bin/env python
"""
Classifies: Dihydroflavonols
Definition: Any hydroxyflavanone in which a hydroxy group is present at position 3 of the heterocyclic ring.
"""

from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol must contain a flavanone (chroman-4-one) core with a free hydroxyl (–OH)
    group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dihydroflavonol, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the flavanone core with atom-mapping:
    #  • The pattern represents a chroman-4-one ring.
    #  • [C;!a:2] and [C;!a:3] are the two saturated carbons.
    #  • We expect that the carbon with map number 3 bears an –OH group.
    core_smarts = "O=C1[C;!a:2][C;!a:3]Oc2ccccc12"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return False, "Internal error: unable to create SMARTS pattern"

    # Find substructure matches for the flavanone core.
    matches = mol.GetSubstructMatches(core_pattern, useChirality=True)
    if not matches:
        return False, "Flavanone (chroman-4-one) core not found"

    # Define a helper to check if an oxygen neighbor looks like a free hydroxyl:
    def has_free_hydroxyl(candidate_atom):
        # For each neighbor of the carbon, if the neighbor is oxygen and it is connected only
        # to the candidate carbon (i.e. not part of the ring linkage which usually connects two heavy atoms),
        # we consider that as a hydroxyl.
        for neighbor in candidate_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                # Check if the bond is a single bond.
                bond = mol.GetBondBetweenAtoms(candidate_atom.GetIdx(), neighbor.GetIdx())
                if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                # A hydroxyl oxygen usually is attached to only one heavy atom.
                # (Note: implicit hydrogens are not counted in GetDegree)
                if neighbor.GetDegree() == 1:
                    return True
        return False

    # Check each occurrence of the core in the molecule
    for match in matches:
        # The match is a tuple of atom indices corresponding to the atoms in core_pattern.
        # We need to extract the atom that corresponds to map 3.
        # To do this, we iterate over the query atoms in core_pattern in order.
        mapped_idx = None
        for i, query_atom in enumerate(core_pattern.GetAtoms()):
            if query_atom.HasProp("molAtomMapNumber"):
                try:
                    map_num = int(query_atom.GetProp("molAtomMapNumber"))
                except ValueError:
                    continue
                if map_num == 3:
                    mapped_idx = match[i]
                    break
        if mapped_idx is None:
            continue  # Should not occur, but proceed safely

        # Now check if the candidate atom (position 3) has a free hydroxyl group.
        candidate_atom = mol.GetAtomWithIdx(mapped_idx)
        if has_free_hydroxyl(candidate_atom):
            return True, "Found flavanone core with a free hydroxyl group at position 3 (dihydroflavonol)."

    return False, "Flavanone core present but no free hydroxyl group found at position 3 of the heterocyclic ring."

# For debugging or unit tests you may add lines like:
# if __name__ == "__main__":
#     test_smiles = "OC1C(Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1"  # Example: (-)-taxifolin
#     result, reason = is_dihydroflavonols(test_smiles)
#     print(result, reason)