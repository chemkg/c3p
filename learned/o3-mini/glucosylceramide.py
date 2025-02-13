"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide 
Definition: Any of the cerebrosides in which the monosaccharide head group is glucose.
This function is_glucosylceramide returns True (with a success message)
if the molecule contains a ceramide backbone with a –CH2OH group that is connected
to a beta-D-glucosyl sugar moiety. Otherwise, it returns False along with a reason.
"""

from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide is a cerebroside having a ceramide backbone with an amide bond
    and a sphingoid base that bears a primary hydroxyl group, to which a beta-D-glucosyl
    sugar is attached.
    
    Our approach is:
     1. Look for a ceramide anchor using a refined SMARTS pattern:
        We require an acyl carbonyl bound to an amide nitrogen, which is in turn
        bound to a chiral carbon that bears a CH2O substructure.
     2. Look for a beta-D-glucosyl moiety using a SMARTS pattern corresponding
        to the pyranose ring of glucose attached via an oxygen (we remove the strict chirality).
     3. Check that one of the oxygen atoms from the sugar is directly bound to the CH2
        (the “attachment” carbon) of the ceramide match.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a glucosylceramide, False otherwise.
      str: Reason for the classification decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a refined SMARTS pattern for the ceramide backbone (anchor):
    # [CX3](=O) matches the acyl carbonyl, then N, then a chiral carbon with a CH2O group.
    ceramide_pattern = Chem.MolFromSmarts("[CX3](=O)N[C@@H](CO)")
    if not ceramide_pattern:
        return False, "Invalid ceramide SMARTS pattern"
    
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if not ceramide_matches:
        return False, "No ceramide backbone (amide motif with CH2O) found"
    
    # Define a SMARTS pattern for a pyranose sugar moiety.
    # Here we look for a six-membered oxygen-containing ring with a pendant hydroxymethyl group.
    # We use a pattern that is not overly strict on chirality.
    glucosyl_pattern = Chem.MolFromSmarts("O1OC(CO)C(O)C(O)C1O")
    if not glucosyl_pattern:
        return False, "Invalid glucosyl SMARTS pattern"
    
    sugar_matches = mol.GetSubstructMatches(glucosyl_pattern)
    if not sugar_matches:
        return False, "No beta-D-glucosyl head group found"
    
    # For a valid glucosylceramide we now require that one of the sugar groups is directly attached
    # to the ceramide backbone. We assume that in our ceramide smart match the fourth atom (index 3)
    # is the CH2 group (which would normally be HO-CH2-) that connects to the sugar.
    for cer_match in ceramide_matches:
        # Get the index of the CH2 carbon in the ceramide match.
        attachment_idx = cer_match[3]
        attachment_atom = mol.GetAtomWithIdx(attachment_idx)
        # Look at all atoms bonded to this CH2 carbon.
        neighbors = [atom.GetIdx() for atom in attachment_atom.GetNeighbors()]
        # Now, check for each sugar match if its first atom (the oxygen expected to link)
        # is directly bonded to the ceramide "attachment" carbon.
        for sugar_match in sugar_matches:
            sugar_attachment_idx = sugar_match[0]
            if sugar_attachment_idx in neighbors:
                return True, "Contains ceramide backbone and beta-D-glucosyl head group (properly connected)"
    
    return False, "Ceramide backbone and beta-D-glucosyl head group found but not properly connected"

# Example usage (you may remove or comment these out in production):
if __name__ == "__main__":
    # Example SMILES for N-octadecanoyl-1-O-beta-D-glucosyl-4-hydroxy-15-methylhexadecasphinganine:
    test_smiles = "CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C"
    result, reason = is_glucosylceramide(test_smiles)
    print(result, reason)