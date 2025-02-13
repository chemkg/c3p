"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: Dihydroflavonols
Definition: Any hydroxyflavanone in which a hydroxy group is present at position 3
of the heterocyclic (chroman-4-one) ring and where position 2 is substituted by an aromatic ring.
The algorithm below first looks for a minimal (non‐chiral) dihydroflavonol core using a SMARTS pattern:
    O[C:2]1[C:3](O)C(=O)C1
This pattern represents the key “chroman-4-one” substructure with:
  • an oxygen (of the pyran ring),
  • a saturated carbon “C:2” (which in the full structure should bear an external aryl substituent),
  • a saturated carbon “C:3” carrying a free –OH (position 3),
  • followed by a carbonyl group.
After a substructure match is found, we check that the atom with mapping number 2 (position 2) has
an external aromatic substituent (i.e. at least one neighbor not part of the match that is aromatic).
"""

from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    
    The algorithm first tries to identify a minimal dihydroflavonol core,
    defined here by the substructure SMARTS "O[C:2]1[C:3](O)C(=O)C1". This pattern
    picks out a chroman-4-one type scaffold in which the saturated carbon (mapped as atom 2)
    is intended to bear an aromatic substituent (i.e. the B-ring). After finding a match,
    we check that the atom with mapping number 2 indeed has at least one neighbor outside the core
    that is aromatic.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a dihydroflavonol, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that captures a minimal dihydroflavonol core.
    # The pattern (ignoring stereochemistry) is meant to represent:
    #   O[C:2]1[C:3](O)C(=O)C1
    # where:
    #   - The first atom (O) is the heterocyclic oxygen.
    #   - Atom with map number 2 (first carbon in ring) is expected to have an aromatic substituent (B-ring).
    #   - Atom with map number 3 is required to be substituted with a free hydroxyl (-OH).
    core_smarts = "O[C:2]1[C:3](O)C(=O)C1"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return False, "Internal error: could not build SMARTS pattern"
    
    # Try to find a substructure match (ignoring chirality)
    matches = mol.GetSubstructMatches(core_pattern, useChirality=False)
    if not matches:
        return False, "Core scaffold (chroman-4-one with 2-aryl substitution and free OH at C3) not found."
    
    # We now check that in at least one match, the atom corresponding to map number 2 has an aryl substituent.
    # To do this we need to know which query atom (by index) in core_pattern was given mapping number "2".
    mapping_index_for_C2 = None
    for i, atom in enumerate(core_pattern.GetAtoms()):
        if atom.HasProp("molAtomMapNumber"):
            if atom.GetProp("molAtomMapNumber") == "2":
                mapping_index_for_C2 = i
                break
    if mapping_index_for_C2 is None:
        return False, "Internal error: mapping for atom C2 not found in pattern."
    
    # Check each match to see if the candidate atom (from mol) corresponding to mapping "2" has an aromatic neighbor.
    for match in matches:
        # match is a tuple; the order corresponds to the query molecule atoms (in the order of core_pattern atoms)
        atom_C2_index = match[mapping_index_for_C2]
        atom_C2 = mol.GetAtomWithIdx(atom_C2_index)
        
        # Look for an external neighbor (not in the core match) that is aromatic.
        has_aromatic_substituent = False
        for nb in atom_C2.GetNeighbors():
            # If the neighbor atom is not part of the core match, check if it is aromatic.
            if nb.GetIdx() not in match and nb.GetIsAromatic():
                has_aromatic_substituent = True
                break
        if has_aromatic_substituent:
            return True, ("Found dihydroflavonol core: chroman-4-one scaffold with an aromatic substituent "
                          "at position 2 and a free hydroxyl at position 3.")
    
    # If no match was found that has an aromatic substituent at the expected position:
    return False, ("Core scaffold (chroman-4-one with 2-aryl substitution and free OH at C3) was found by SMARTS "
                   "but the candidate atom at position 2 does not appear to be substituted by an aromatic ring.")

# Uncomment the following for local testing:
# if __name__ == "__main__":
#     test_smiles = "O[C@H]1[C@@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1"  # Example: (-)-taxifolin
#     result, reason = is_dihydroflavonols(test_smiles)
#     print(result, reason)