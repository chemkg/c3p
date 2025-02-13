"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide
Definition: Any of the cerebrosides in which the monosaccharide head group is glucose.
This function is_glucosylceramide returns True if the molecule contains a ceramide backbone 
(with an N-acyl sphingoid fragment featuring a CH2 group that connects via an oxygen bridge to 
a six‐membered ring having one oxygen and five carbons, as expected for a glucopyranose sugar).
"""

from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    Our approach is:
      1. Look for a ceramide backbone defined by an acyl group attached to an amide nitrogen
         and then to a CH (which bears a CH2O fragment). The SMARTS pattern we use here is:
         "[CX3](=O)N[CH](CO)".
      2. From the attachment atom (the CH2 that is part of the "CO" fragment), check its
         oxygen neighbors. We require that at least one such oxygen is linked (directly or via
         a short bridge) to a six-membered ring.
      3. Inspect that ring to verify it is composed of exactly one oxygen and five carbons,
         consistent with a glucopyranose ring.
      
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a glucosylceramide, False otherwise.
      str: Explanation for the classification decision.
    """
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern to match a typical ceramide fragment.
    # This looks for: a carbonyl carbon (acyl function) attached to an N which is then attached
    # to a carbon (which should carry a –CH2O fragment as further elaborated in the SMILES).
    ceramide_pattern = Chem.MolFromSmarts("[CX3](=O)N[CH](CO)")
    if ceramide_pattern is None:
        return False, "Invalid ceramide SMARTS pattern"
    
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if not ceramide_matches:
        return False, "No ceramide backbone (C(=O)N[CH](CO)) found"
    
    # Get ring information from the molecule for sugar ring evaluation.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices.
    
    # Loop over each found ceramide match.
    # In our SMARTS the match is expected to return 4 atoms:
    #   index 0: acyl carbonyl carbon,
    #   index 1: amide N,
    #   index 2: the chiral carbon (sphingoid base part),
    #   index 3: the carbon in the "CO" fragment (attachment point to the sugar)
    for match in ceramide_matches:
        if len(match) < 4:
            continue  # Skip matches that are not full.
        attachment_idx = match[3]  # The attachment carbon (CH2 fragment).
        attachment_atom = mol.GetAtomWithIdx(attachment_idx)
        # Look at oxygen neighbors of the attachment atom.
        for neighbor in attachment_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 8:  # Only consider oxygen atoms.
                continue
            oxy_idx = neighbor.GetIdx()
            # Now check if this oxygen is part of a six-membered ring.
            for ring in rings:
                if oxy_idx in ring and len(ring) == 6:
                    # Count the atoms in the ring.
                    o_count = 0
                    c_count = 0
                    for idx in ring:
                        atom = mol.GetAtomWithIdx(idx)
                        if atom.GetAtomicNum() == 8:
                            o_count += 1
                        elif atom.GetAtomicNum() == 6:
                            c_count += 1
                    # A glucopyranose ring should have exactly 1 oxygen and 5 carbons.
                    if o_count == 1 and c_count == 5:
                        return True, ("Found a ceramide backbone (with a C(=O)N[CH](CO) fragment) " +
                                      "whose CH2O attachment is linked via an oxygen to a 6-membered " +
                                      "ring (1 oxygen, 5 carbons); indicative of a beta-D-glucosyl moiety.")
                        
    return False, ("Ceramide backbone and sugar attachment found, but no appropriate " +
                   "beta-D-glucosyl sugar ring (6-membered ring with 1 oxygen and 5 carbons) detected")

# Example usage:
if __name__ == "__main__":
    # Test the function with one of the provided glucosylceramide SMILES:
    test_smiles = ("CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)"
                   "[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C")
    result, reason = is_glucosylceramide(test_smiles)
    print(result, reason)