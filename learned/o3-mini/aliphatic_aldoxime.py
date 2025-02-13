"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: Any aldoxime derived from an aliphatic aldehyde.
"""

from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is any aldoxime (C=N–OH) that is derived from an aliphatic aldehyde.
    This means that the carbon of the aldoxime (originally the aldehyde carbon) should be a CH unit,
    and should not be directly or indirectly connected to an aromatic system.
    
    Our strategy is to identify the aldoxime group using a SMARTS query and then
    verify that the carbon involved (1) is not aromatic and (2) none of its nearby heavy-atom neighbors
    (apart from the bonded nitrogen) are either aromatic themselves or directly attached to aromatic atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an aliphatic aldoxime, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so we can check for CH1
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for the aldoxime moiety:
    # [CH1;!a]  ensures a CH unit that is not flagged aromatic, attached by double bond to [N],
    # which in turn is single bonded to [O;H] (an oxygen with an explicit hydrogen).
    aldoxime_smarts = "[CH1;!a]=[N][O;H]"
    aldoxime_pattern = Chem.MolFromSmarts(aldoxime_smarts)
    if aldoxime_pattern is None:
        return False, "Error in defining SMARTS for aldoxime"
    
    # Find matches for the aldoxime group.
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "No aldoxime group (C=N–OH with a CH unit) found"
    
    # For each match, perform extra aliphatic filters:
    for match in matches:
        # The SMARTS pattern yields indices for (carbon, nitrogen, oxygen).
        carbon_idx, nitrogen_idx, oxygen_idx = match
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # First check: The carbon atom should not be aromatic.
        if carbon_atom.GetIsAromatic():
            # Not aliphatic if the key carbon itself is aromatic.
            continue
        
        # Now check neighbors of the aldehyde carbon (excluding the nitrogen that forms the oxime).
        qualifies = True
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIdx() == nitrogen_idx:
                continue  # skip the nitrogen of the oxime group
            if neighbor.GetAtomicNum() == 1:
                continue  # skip hydrogens
            
            # Direct neighbor check: if neighbor itself is aromatic, reject.
            if neighbor.GetIsAromatic():
                qualifies = False
                break
            
            # Second-level check: look at the neighbors of this neighbor (excluding the original carbon).
            # If any neighbor (of the neighbor) is aromatic, then the substituent is benzylic to an aromatic ring.
            for sub_neighbor in neighbor.GetNeighbors():
                if sub_neighbor.GetIdx() == carbon_atom.GetIdx():
                    continue
                if sub_neighbor.GetAtomicNum() != 1 and sub_neighbor.GetIsAromatic():
                    qualifies = False
                    break
            if not qualifies:
                break
        
        if qualifies:
            return True, "Aldoxime group derived from an aliphatic aldehyde found"
    
    return False, "Aldoxime group found but appears to be derived from an aromatic aldehyde"

# Optional: For testing purposes, you can call the function on one example.
if __name__ == "__main__":
    test_smiles = "[H]\\C(C)=N\\O"  # (Z)-acetaldehyde oxime
    result, reason = is_aliphatic_aldoxime(test_smiles)
    print(result, "->", reason)