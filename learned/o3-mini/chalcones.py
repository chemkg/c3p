"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: chalcones (1,3-diphenylpropenone and its substituted derivatives)
A chalcone is defined as a ketone with the core α,β‐unsaturated system attached
to two aromatic rings (Ar–CH=CH–C(=O)–Ar or its reverse).
The program checks for one of these substructures and further verifies that the
unsaturated “bridge” (CH=CH and C(=O) atoms) is exocyclic (non-ring) while the flanking
aromatic atoms are in rings.
"""

from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone (or derivative) based on its SMILES string.
    A chalcone contains an alpha,beta-unsaturated ketone linker connecting two aromatic rings.
    
    This function uses two SMARTS patterns corresponding to the two possible connectivities:
      Pattern 1: Aromatic–CH=CH–C(=O)–Aromatic
      Pattern 2: Aromatic–C(=O)–CH=CH–Aromatic

    It then further checks that in any match:
      • The flanking atoms (first and last in the match) are aromatic and are in a ring.
      • The bridging atoms (the two alkene carbons and the carbonyl carbon) are not in a ring.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule contains a chalcone-like substructure, False otherwise.
        str: A descriptive reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the chalcone motifs.
    # Pattern 1: aromatic–CH=CH–C(=O)–aromatic
    chalcone_pattern1 = Chem.MolFromSmarts("a-[CH]=[CH]-[C](=O)-a")
    # Pattern 2: aromatic–C(=O)–CH=CH–aromatic (reverse connectivity)
    chalcone_pattern2 = Chem.MolFromSmarts("a-[C](=O)-[CH]=[CH]-a")
    
    # Helper function: Given a substructure match (list of atom indices),
    # check that the endpoints are aromatic and in a ring and the bridging atoms are exocyclic.
    def valid_chalcone_match(match):
        # We expect match to have 5 atoms:
        # For both patterns:
        #   match[0]: aromatic atom (should be in ring)
        #   match[1]: first bridging atom (should not be in any ring)
        #   match[2]: second bridging atom (should not be in any ring)
        #   match[3]: carbonyl carbon (should not be in any ring)
        #   match[4]: aromatic atom (should be in ring)
        if len(match) != 5:
            return False
        # Check endpoints: must be aromatic and in a ring.
        for idx in [match[0], match[4]]:
            atom = mol.GetAtomWithIdx(idx)
            if not (atom.GetIsAromatic() and atom.IsInRing()):
                return False
        # Check bridging atoms: should not be in any ring.
        for idx in match[1:4]:
            if mol.GetAtomWithIdx(idx).IsInRing():
                return False
        return True
    
    found = False
    reason_details = ""
    
    # Look for pattern 1 matches.
    matches1 = mol.GetSubstructMatches(chalcone_pattern1)
    for match in matches1:
        if valid_chalcone_match(match):
            found = True
            break  # Once a valid match is found we can stop checking.
    
    # If not found, try pattern 2.
    if not found:
        matches2 = mol.GetSubstructMatches(chalcone_pattern2)
        for match in matches2:
            if valid_chalcone_match(match):
                found = True
                break
    
    if found:
        return True, "Molecule contains a chalcone-like α,β-unsaturated ketone substructure with two aromatic rings."
    else:
        return False, "Molecule does not contain the chalcone core (α,β-unsaturated ketone linked to two aromatic rings) as defined."

# Example usage for testing:
if __name__ == "__main__":
    # Testing with cis-chalcone as example.
    smiles_example = "O=C(c1ccccc1)\\C=C/c1ccccc1"
    result, reason = is_chalcones(smiles_example)
    print(f"Result: {result}\nReason: {reason}")