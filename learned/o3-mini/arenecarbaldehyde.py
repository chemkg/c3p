"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: Arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is directly attached to an arene.
An arenecarbaldehyde must have a carbonyl group (C=O) where the carbon is bonded to an aromatic carbon,
with that aromatic atom being part of a ring that contains only carbon atoms.
"""

from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    
    A valid arenecarbaldehyde contains a carbonyl group typical of an aldehyde where:
      - The carbonyl carbon is sp2 and bonded to an oxygen (as C=O) and only one other heavy atom
        (the aromatic neighbor); the remaining valence comes from an implicit hydrogen.
      - The aromatic neighbor is a carbon (atomic number 6) that is aromatic and
        belongs to at least one aromatic ring made up entirely of carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is an arenecarbaldehyde, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule is aromatized properly.
    Chem.SanitizeMol(mol)
    
    # SMARTS pattern to match a carbonyl attached to an aromatic carbon.
    # Note: We have dropped the explicit hydrogen requirement so that implicit hydrogens are allowed.
    aldehyde_arene_smarts = "[CX3](=O)[c]"
    pattern = Chem.MolFromSmarts(aldehyde_arene_smarts)
    if pattern is None:
        return False, "Error creating SMARTS pattern"
    
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No aldehyde group directly attached to an aromatic atom was found"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # For each aldehyde match, check that:
    # 1. The carbonyl carbon behaves like an aldehyde (only one heavy neighbor aside from oxygen).
    # 2. Its aromatic neighbor is indeed an aromatic carbon and belongs to a ring that contains only carbons.
    for match in matches:
        # match[0] is the carbonyl carbon, match[1] is the attached aromatic neighbor.
        carbonyl_c = mol.GetAtomWithIdx(match[0])
        aromatic_neighbor = mol.GetAtomWithIdx(match[1])
        
        # Check that the neighbor is a carbon atom.
        if aromatic_neighbor.GetAtomicNum() != 6:
            continue
        
        # Check that the carbonyl carbon is consistent with an aldehyde:
        # It should have two explicit bonds (to the oxygen and to the aromatic neighbor)
        # and one implicit hydrogen (totalling a degree of 3).
        if carbonyl_c.GetTotalDegree() != 3:
            # If there is more than one heavy atom attached (excluding oxygen), it is likely a ketone.
            continue
        
        # Next, confirm that the aromatic neighbor is marked as aromatic.
        if not aromatic_neighbor.GetIsAromatic():
            continue
        
        # Now check that the aromatic neighbor is part of at least one ring that is a pure carbon arene.
        neighbor_idx = aromatic_neighbor.GetIdx()
        pure_arene_found = False
        for ring in ring_info:
            if neighbor_idx in ring:
                # Check that every atom in the ring that is aromatic is carbon.
                all_carbons = True
                for idx in ring:
                    a = mol.GetAtomWithIdx(idx)
                    if a.GetIsAromatic() and a.GetAtomicNum() != 6:
                        all_carbons = False
                        break
                if all_carbons:
                    pure_arene_found = True
                    break
        if pure_arene_found:
            return True, ("Found an aldehyde group with the carbonyl carbon (with an implicit hydrogen) " 
                          "attached to an aromatic carbon that belongs to a pure arene ring (all aromatic atoms are carbon).")
    
    return False, "Aldehyde group found, but the attached aromatic atom is not part of a pure arene ring or does not show correct aldehyde bonding."
    
# Example usage (uncomment to test):
# test_smiles = [
#     "[H]C(=O)c1cccc2ccccc12",  # 1-naphthaldehyde (should be True)
#     "OC(=O)c1cc(O)c(C(=O)c2c(O)cccc2O)c(C=O)c1",  # TAN-931 (should be True)
#     "Cc1cccc(C=O)c1",          # m-tolualdehyde (should be True)
#     "[O-][N+](=O)c1ccc(C=O)cc1",  # 4-nitrobenzaldehyde (should be True)
# ]
# for smi in test_smiles:
#     result, reason = is_arenecarbaldehyde(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")