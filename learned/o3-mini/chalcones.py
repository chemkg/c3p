"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: chalcones 
A chalcone is defined as a ketone with a 1,3-diphenylpropenone core,
i.e. an α,β‐unsaturated ketone with two aromatic (ring) substituents. The core can have either connectivity:
  Pattern 1: Aromatic–CH=CH–C(=O)–Aromatic
  Pattern 2: Aromatic–C(=O)–CH=CH–Aromatic
This program uses two SMARTS queries and then verifies that:
  • The two outer atoms are aromatic and part of a ring.
  • The bridging atoms (vinyl and carbonyl carbons) are exocyclic (not in a ring).
  • The expected carbonyl carbon shows a double bond to an oxygen.
"""

from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone (or derivative) based on its SMILES string.
    A chalcone contains an α,β‐unsaturated ketone linker connecting two aromatic rings.
    
    The function uses two SMARTS patterns:
      Pattern 1: Aromatic–CH=CH–C(=O)–Aromatic
      Pattern 2: Aromatic–C(=O)–CH=CH–Aromatic
      
    For each match we then check:
      • The endpoints (first and last atoms) are aromatic and in a ring.
      • The bridging atoms (the vinyl/carbonyl carbons) are not in any ring.
      • The carbonyl carbon shows a double bond to an oxygen.
      • Also, the expected C=C double bond is present in the appropriate position.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if a valid chalcone core is found, False otherwise.
        str: A descriptive reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for the two possible connectivities.
    # Pattern 1: Aromatic–C=C–C(=O)–Aromatic
    pattern1 = Chem.MolFromSmarts("a-[C]=[C]-[C](=O)-a")
    # Pattern 2: Aromatic–C(=O)–C=C–Aromatic
    pattern2 = Chem.MolFromSmarts("a-[C](=O)-[C]=[C]-a")
    
    # Helper to check if an atom (by index) has a double bond to an oxygen.
    def has_double_bonded_oxygen(atom_idx):
        atom = mol.GetAtomWithIdx(atom_idx)
        for bond in atom.GetBonds():
            # Check if bond is double and the neighboring atom is oxygen.
            if bond.GetBondTypeAsDouble() == 2:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    return True
        return False

    # Helper to validate a single match.
    # match is a tuple/list of 5 atom indices.
    def valid_match(match):
        if len(match) != 5:
            return False
        
        # Check endpoints: must be aromatic and in a ring.
        for idx in (match[0], match[4]):
            atom = mol.GetAtomWithIdx(idx)
            if not (atom.GetIsAromatic() and atom.IsInRing()):
                return False
        
        # Check bridging atoms (indices 1,2,3) are exocyclic (not in a ring).
        for idx in (match[1], match[2], match[3]):
            if mol.GetAtomWithIdx(idx).IsInRing():
                return False
        
        # Now determine if the match fits pattern1 or pattern2 by checking for the correct double bonds.
        # For pattern1 [a]-[C]=[C]-[C](=O)-[a]:
        #   => the bond between atoms at match[1] and match[2] must be double,
        #      and the carbonyl (match[3]) must have a double bond to oxygen.
        bond_12 = mol.GetBondBetweenAtoms(match[1], match[2])
        pattern1_valid = (bond_12 is not None and bond_12.GetBondTypeAsDouble() == 2 and
                          has_double_bonded_oxygen(match[3]))
        # For pattern2 [a]-[C](=O)-[C]=[C]-[a]:
        #   => the bond between atoms at match[2] and match[3] must be double,
        #      and the carbonyl (match[1]) must have a double bond to oxygen.
        bond_23 = mol.GetBondBetweenAtoms(match[2], match[3])
        pattern2_valid = (bond_23 is not None and bond_23.GetBondTypeAsDouble() == 2 and
                          has_double_bonded_oxygen(match[1]))
        
        if pattern1_valid or pattern2_valid:
            return True
        return False

    # Try matching pattern1.
    matches = mol.GetSubstructMatches(pattern1)
    for match in matches:
        if valid_match(match):
            return True, "Molecule contains a chalcone-like α,β-unsaturated ketone substructure (Pattern 1)."
    
    # If nothing found yet, try pattern2.
    matches = mol.GetSubstructMatches(pattern2)
    for match in matches:
        if valid_match(match):
            return True, "Molecule contains a chalcone-like α,β-unsaturated ketone substructure (Pattern 2)."
    
    return False, "Molecule does not contain the chalcone core (α,β-unsaturated ketone linked to two aromatic rings) as defined."

# Example usage for testing:
if __name__ == "__main__":
    # Example: cis-chalcone should return True.
    test_smiles = "O=C(c1ccccc1)\\C=C/c1ccccc1"
    result, reason = is_chalcones(test_smiles)
    print(f"Test SMILES: {test_smiles}\nResult: {result}\nReason: {reason}")