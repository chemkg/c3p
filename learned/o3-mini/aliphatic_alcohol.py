"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: Aliphatic Alcohol (an alcohol derived from an aliphatic compound)

An aliphatic alcohol is defined here as one containing an –OH group that is
attached to a non‐aromatic (sp3) carbon. In addition to using a SMARTS pattern,
we further require that the sp3 carbon bearing the –OH is not directly bonded
to any aromatic carbons. This extra check helps avoid false positives for molecules
that contain an alcohol group in an otherwise complex (e.g. glycosidic or aromatic)
environment.
"""

from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    
    An aliphatic alcohol is defined as one that contains at least one hydroxyl (-OH)
    group attached to a non-aromatic (sp3) carbon. In addition, we require that the
    sp3 carbon (directly bound to the -OH) is not directly attached to an aromatic carbon,
    so as to avoid false positives from complex glycosides or molecules that contain
    aromatic rings.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an aliphatic alcohol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for an alcohol where -OH is attached
    # to a sp3 aliphatic carbon. Here we use [#6X4;!$([#6](=O))] to indicate a
    # saturated carbon that is not part of a carbonyl group, followed by [OX2H].
    alcohol_smarts = "[#6X4;!$([#6](=O))][OX2H]"
    aliphatic_alcohol_pattern = Chem.MolFromSmarts(alcohol_smarts)
    
    # First, find all substructure matches for the alcohol pattern.
    matches = mol.GetSubstructMatches(aliphatic_alcohol_pattern)
    if not matches:
        return False, "No -OH group attached to a sp3 (non-aromatic) carbon found"

    # Now check the neighborhood of each match.
    # Each match is a tuple of atom indices corresponding to the SMARTS pattern.
    # Convention: match[0] is the carbon atom, match[1] is the O (with -H).
    valid_match_found = False
    for match in matches:
        carbon_idx, oxygen_idx = match[0], match[1]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        # If the carbon is marked as aromatic (should not be since we used #6X4),
        # then skip this match.
        if carbon_atom.GetIsAromatic():
            continue
        # Check the immediate neighbors of the carbon (except the -OH oxygen)
        adjacent_aromatic = False
        for neighbor in carbon_atom.GetNeighbors():
            # Skip the oxygen that is part of the alcohol group
            if neighbor.GetIdx() == oxygen_idx:
                continue
            if neighbor.GetIsAromatic():
                adjacent_aromatic = True
                break
        if adjacent_aromatic:
            continue
        # If we get here then this alcohol group is truly attached to an aliphatic environment.
        valid_match_found = True
        break
    
    if not valid_match_found:
        return False, ("Found an alcohol group attached to a sp3 carbon, "
                       "but the carbon is directly bonded to an aromatic system")
    
    # If at least one match passes, we classify the molecule as an aliphatic alcohol.
    return True, "Molecule contains an aliphatic alcohol group: -OH is attached to a non-aromatic sp3 carbon"

# Example usage: (These can be commented out when used as a module)
if __name__ == '__main__':
    example_smiles = [
        "OC(C(C)C)CC/C=C/C=C/C([C@@H](O)C)(C)C",  # Graphostromol G (valid)
        "CCCCCCCCCCCCC(O)CCCC",                  # octadecan-5-ol (valid)
        "CO",                                    # methanol (valid)
        "C[C@@H](C1=CC=CC=C1)NC(=O)[C@@H]2[C@H](...)",  # false positive example (has aromatic ring)
        "c1ccccc1O"                             # phenol (invalid: -OH on aromatic carbon)
    ]
    for smi in example_smiles:
        result, reason = is_aliphatic_alcohol(smi)
        print(f"SMILES: {smi}\n  Classified as aliphatic alcohol? {result}\n  Reason: {reason}\n")