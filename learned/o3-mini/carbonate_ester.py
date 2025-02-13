"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: Carbonate Ester
Definition: Any carbonate that is carbonic acid in which the hydrogens have been replaced by organyl groups.
That is, a carbonate ester is a derivative of carbonic acid (O=C(OH)2) in which at least one -OH group has been replaced
by an organo (typically –OC(alkyl)) group. This includes diesters (e.g. dimethyl carbonate), monoesters (e.g. monoethyl carbonate)
as well as cyclic carbonate derivatives (e.g. ethylene carbonate).
"""

from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    
    The method looks for a central carbon atom that is connected in a carbonate fashion:
      - It is double bonded to one oxygen.
      - It is single bonded to two oxygens.
      - In the parent carbonic acid the two singly-bound oxygens would have an H each.
      - In a carbonate ester, at least one of these oxygens is substituted with a non-H (an organyl group).
    
    To enforce that we are dealing with an unmodified carbonic acid scaffold (and not a more extended acyl unit),
    we require that the central carbon (the one bearing the carbonate motif) is bonded to exactly three heavy atoms.
    
    We handle three cases via SMARTS:
      (a) A diester pattern: both singly-bound oxygens are substituted by a carbon (or, more generally, a non-H).
      (b) Two monoester patterns: one version where the first oxygen is unsubstituted (i.e. remains –OH)
          and the second is substituted and vice versa.
      (c) Cyclic carbonates are expected to satisfy the diester condition because in the ring the oxygen(s) are bound to carbons.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens so that -OH groups become visible
    mol = Chem.AddHs(mol)
    
    # Define SMARTS patterns with explicit atom mappings.
    # The central carbonate carbon is tagged as :1.
    # For a diester, both oxygens (tagged :2 and :3) are bound to a heavy atom (here we require a bound carbon, [#6])
    patt_diester = Chem.MolFromSmarts("[C:1](=O)([O:2][#6])([O:3][#6])")
    
    # For a monoester, one of the singly-bound oxygens is an –OH (explicit H attached, so [OX2H])
    patt_monoester1 = Chem.MolFromSmarts("[C:1](=O)([O:2][H])([O:3][#6])")
    patt_monoester2 = Chem.MolFromSmarts("[C:1](=O)([O:2][#6])([O:3][H])")
    
    # List of SMARTS candidate patterns
    patterns = [
        ("diester", patt_diester),
        ("monoester", patt_monoester1),
        ("monoester", patt_monoester2)
    ]
    
    # For each pattern, see if it appears in the molecule.
    for patt_type, patt in patterns:
        matches = mol.GetSubstructMatches(patt)
        if matches:
            # For every match, do an extra check on the central carbon.
            for match in matches:
                # match is a tuple of atom indices corresponding to the mapped atoms [C:1], [O:2], [O:3]
                carbon_idx = match[0]
                carbon_atom = mol.GetAtomWithIdx(carbon_idx)
                # Count heavy-atom (non-hydrogen) neighbors of the central carbon.
                heavy_neighbors = [nbr for nbr in carbon_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if len(heavy_neighbors) != 3:
                    # This likely means the carbon is further substituted; skip this match.
                    continue
                # If we get here, we have a good match for a carbonate ester.
                if patt_type == "diester":
                    reason = "Found carbonate moiety (O=C(–O–R)2) with both –OH replaced by organyl groups."
                else:
                    reason = "Found carbonate moiety (O=C(–OH)(–O–R)) with at least one organyl substitution."
                return True, reason
    # If no match was found, either no carbonate motif is present or it did not satisfy the restraints.
    return False, "No carbonate ester motif (O=C(–OH/–O–R)2 with correct connectivity) detected in the molecule."

# For quick testing, you can uncomment the lines below.
# test_smiles = [
#     "O(CCC(C)=C)C(OCC)=O",  # Ethyl 3-methylbut-3-enyl carbonate
#     "COC(=O)OC",            # dimethyl carbonate
#     "O1CCOC1=O",            # Ethylene carbonate
#     "CCOC(O)=O",            # monoethyl carbonate
#     "O=C(O)O",              # carbonic acid (should be False)
# ]
# for smi in test_smiles:
#     result, reason = is_carbonate_ester(smi)
#     print(f"SMILES: {smi} --> {result}. Reason: {reason}")