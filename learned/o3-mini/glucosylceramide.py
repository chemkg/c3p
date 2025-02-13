"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide
Definition: Any of the cerebrosides in which the monosaccharide head group is glucose.
This function is_glucosylceramide returns True if the molecule contains a ceramide backbone 
(with an amide group and a CH2O “attachment” site) that is linked via an oxygen to a six‐membered 
ring consisting of one oxygen and five carbons (as in a glucopyranose ring).
"""

from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    Our approach is:
      1. Identify a ceramide backbone by looking for an acyl group 
         (carbonyl [CX3](=O)) connected to an amide nitrogen and then to a –CH2 group.
      2. From that CH2 group (our presumed sugar attachment point), check its neighboring atoms.
         Look for an oxygen that is part of a six‐membered ring.
      3. In that ring, verify that exactly one oxygen is present (the typical composition 
         of a glucopyranose ring) and five carbons.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a glucosylceramide, False otherwise.
      str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a ceramide backbone.
    # We look for an acyl carbonyl [CX3](=O) attached to an N and then to a carbon bearing a CH2(O...) fragment.
    # We are not overly strict with chiral specifications to accommodate variations.
    ceramide_pattern = Chem.MolFromSmarts("[CX3](=O)N[C](CO)")
    if ceramide_pattern is None:
        return False, "Invalid ceramide SMARTS pattern"
    
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if not ceramide_matches:
        return False, "No ceramide backbone (amide bonded to a CH2O group) found"
    
    # Get ring information for later sugar ring evaluation.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # tuple of atom index tuples
    
    # Loop over each matched ceramide backbone
    for match in ceramide_matches:
        # We assume our pattern matches four atoms: 
        #   0: carbonyl carbon, 1: carbonyl oxygen is implicit in pattern, 1: N, 2: next carbon, 3: the CH2 (attachment) group.
        if len(match) < 4:
            continue
        attachment_idx = match[3]
        attachment_atom = mol.GetAtomWithIdx(attachment_idx)
        # Check neighbors of the attachment carbon to see if it is linked to an oxygen
        for neighbor in attachment_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 8:  # only consider oxygen atoms
                continue
            oxy_idx = neighbor.GetIdx()
            # Look for a ring containing this oxygen
            for ring in rings:
                if oxy_idx in ring and len(ring) == 6:
                    # Count how many oxygens are in this ring
                    o_count = 0
                    c_count = 0
                    for idx in ring:
                        atom = mol.GetAtomWithIdx(idx)
                        if atom.GetAtomicNum() == 8:
                            o_count += 1
                        elif atom.GetAtomicNum() == 6:
                            c_count += 1
                    # A typical glucopyranose ring should have 1 oxygen and 5 carbons.
                    if o_count == 1 and c_count == 5:
                        return True, ("Found a ceramide backbone with a CH2O attachment "
                                      "linked to a 6-membered (pyranose) sugar ring (indicative of beta-D-glucosyl).")
    return False, "Ceramide backbone and sugar attachment found, but no appropriate beta-D-glucosyl sugar ring detected"

# Example usage (you may remove or comment these out in production):
if __name__ == "__main__":
    # Example SMILES for N-octadecanoyl-1-O-beta-D-glucosyl-4-hydroxy-15-methylhexadecasphinganine
    test_smiles = "CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C"
    result, reason = is_glucosylceramide(test_smiles)
    print(result, reason)