"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: chalcones
A chalcone is defined as a ketone with a 1,3-diphenylpropenone core,
i.e. an α,β-unsaturated ketone with two aromatic (ring) substituents.
Accepted connectivities are:
  Pattern 1: Aromatic–CH=CH–C(=O)–Aromatic
  Pattern 2: Aromatic–C(=O)–CH=CH–Aromatic

This program uses two SMARTS queries with improved features:
  • SMARTS now use "[cR]" for aromatic ring carbons ensuring proper matching
    even when substitution is present.
  • The bridging atoms are checked to be exocyclic (not part of any ring).
  • The carbonyl carbon is verified to have a double bond to an oxygen.
"""

from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone (or derivative) based on its SMILES string.
    A chalcone is defined as containing an α,β‐unsaturated ketone linker connecting
    two aromatic rings. Two connectivities are allowed:
      Pattern 1: Aromatic–CH=CH–C(=O)–Aromatic
      Pattern 2: Aromatic–C(=O)–CH=CH–Aromatic
    The SMARTS queries use "[cR]" for aromatic ring carbons.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if a valid chalcone core is found, False otherwise.
        str: A descriptive reason for the classification.
    """
    
    # Parse SMILES, return an error if invalid.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for detecting chalcone cores.
    # Pattern 1: Aromatic–CH=CH–C(=O)–Aromatic
    pattern1 = Chem.MolFromSmarts("[cR]-[#6]=[#6]-[#6](=O)-[cR]")
    # Pattern 2: Aromatic–C(=O)–CH=CH–Aromatic
    pattern2 = Chem.MolFromSmarts("[cR]-[#6](=O)-[#6]=[#6]-[cR]")
    
    # Helper function that checks if an atom (by index) has a double-bonded oxygen.
    def has_double_bonded_oxygen(atom_idx):
        atom = mol.GetAtomWithIdx(atom_idx)
        for bond in atom.GetBonds():
            # Check for a double bond (bond order equal to 2) to an oxygen.
            if bond.GetBondTypeAsDouble() == 2:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    return True
        return False

    # Helper function that verifies if a match (a tuple of atom indices) satisfies our conditions.
    # For both patterns, the match is expected to be a tuple of 5 atoms.
    # In pattern1:
    #    match[0] and match[4] (the endpoints) are expected to be aromatic and in rings,
    #    match[1], match[2], match[3] (the bridging atoms) must be exocyclic.
    #    The carbonyl carbon is at index 3 and must have a double-bonded oxygen.
    # In pattern2:
    #    match[0] and match[4] are endpoints,
    #    match[1], match[2], match[3] are bridging (with match[1] being the carbonyl in this case).
    def valid_match(match, pattern_variant):
        # Ensure we got 5 atoms.
        if len(match) != 5:
            return False
        
        # Check endpoints: they must be aromatic and be part of a ring.
        for idx in (match[0], match[4]):
            atom = mol.GetAtomWithIdx(idx)
            if not (atom.GetIsAromatic() and atom.IsInRing()):
                return False
        
        # Check bridging atoms (should be exocyclic, i.e. not in a ring).
        # In both patterns positions 1,2,3 must not be in any ring.
        for idx in (match[1], match[2], match[3]):
            if mol.GetAtomWithIdx(idx).IsInRing():
                return False
        
        # Check connectivity specific to the variant.
        if pattern_variant == 1:
            # For pattern1: the bond between atoms at indices 1 and 2 must be double.
            bond_12 = mol.GetBondBetweenAtoms(match[1], match[2])
            if bond_12 is None or bond_12.GetBondTypeAsDouble() != 2:
                return False
            # Verify that the carbonyl carbon (match[3]) has a double-bonded oxygen.
            if not has_double_bonded_oxygen(match[3]):
                return False
        elif pattern_variant == 2:
            # For pattern2: the bond between atoms at indices 2 and 3 must be double.
            bond_23 = mol.GetBondBetweenAtoms(match[2], match[3])
            if bond_23 is None or bond_23.GetBondTypeAsDouble() != 2:
                return False
            # Verify that the carbonyl carbon (match[1]) has a double-bonded oxygen.
            if not has_double_bonded_oxygen(match[1]):
                return False
        else:
            return False
        
        return True

    # Try matching with pattern1.
    matches1 = mol.GetSubstructMatches(pattern1)
    for match in matches1:
        if valid_match(match, pattern_variant=1):
            return True, "Molecule contains a chalcone-like α,β-unsaturated ketone substructure (Pattern 1)."
    
    # If pattern1 did not match, try pattern2.
    matches2 = mol.GetSubstructMatches(pattern2)
    for match in matches2:
        if valid_match(match, pattern_variant=2):
            return True, "Molecule contains a chalcone-like α,β-unsaturated ketone substructure (Pattern 2)."
    
    return False, "Molecule does not contain the chalcone core (α,β-unsaturated ketone linked to two aromatic rings) as defined."

# Example usage for testing:
if __name__ == "__main__":
    # Test a canonical chalcone (cis-chalcone).
    test_smiles = "O=C(c1ccccc1)\\C=C/c1ccccc1"
    result, reason = is_chalcones(test_smiles)
    print(f"Test SMILES: {test_smiles}\nResult: {result}\nReason: {reason}")