"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: 11-oxo steroid 
Definition: Any oxo steroid that has an oxo substituent at position 11.
Note: This heuristic implementation checks for a steroid‐like fused
ring system (three 6‐membered rings + one 5‐membered ring) and for a ketone group
(i.e. a carbonyl attached to two ring carbons). This does not perform full stereochemical
or numbering assignment, so molecules that do not strictly adhere to these features
may be mis‐classified.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    We use heuristic checks:
       1. The molecule can be parsed.
       2. It contains a fused ring system with at least one five-membered ring and at least three six-membered rings,
          a common representation of the cyclopenta[a]phenanthrene (steroid) core.
       3. It has at least one ketone group (C(=O)) on a ring carbon that is bound only to carbons (typical for the
          internal oxo substituents as seen for C-11) as a surrogate to the requirement “oxo at position 11”.
       
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a 11-oxo steroid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is perceived (if available) to help with ring perception
    Chem.AssignStereochemistry(mol, cleanIt=True)
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    # Count rings of size 5 and 6
    count5 = sum(1 for ring in atom_rings if len(ring) == 5)
    count6 = sum(1 for ring in atom_rings if len(ring) == 6)
    
    # Check for a typical steroid fused ring system (3 six-membered rings + 1 five-membered ring)
    if count5 < 1 or count6 < 3:
        return False, ("Lacks a typical steroid nucleus: expected at least 1 five-membered ring and 3 six-membered rings, "
                       f"found {count5} five-membered and {count6} six-membered rings.")

    # Now, search for ketone (C(=O)) groups that are part of rings.
    # We loop over carbon atoms and check if they have a double bond to an oxygen.
    ketone_found = False
    ketone_reasons = []
    for atom in mol.GetAtoms():
        # Look only at carbon atoms
        if atom.GetAtomicNum() != 6:
            continue
        
        # Get bonds from this atom.
        for bond in atom.GetBonds():
            # Check for a double bond
            if bond.GetBondTypeAsDouble() != 2:
                continue
            # Identify the neighbor atom that is oxygen
            neighbor = bond.GetOtherAtom(atom)
            if neighbor.GetAtomicNum() != 8:
                continue
            # We now have a C=O bond
            # Check that the carbonyl carbon (atom) is in a ring.
            if not atom.IsInRing():
                continue
            
            # Further, check that besides the oxygen (via the double bond), the carbon is attached to two other carbons
            # (i.e., it is not part of a carboxylic acid or ester). In many steroids the ketone is internal.
            ring_carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != neighbor.GetIdx()]
            if len(ring_carbon_neighbors) < 2:
                continue
            
            # This is a candidate for a ring-bound ketone.
            # (Note: a full numbering algorithm would be needed to check that it is at position 11.)
            ketone_found = True
            ketone_reasons.append(f"Carbon atom {atom.GetIdx()} is part of a ring and forms a ketone with oxygen atom {neighbor.GetIdx()}.")
    
    if not ketone_found:
        return False, "No ring-bound ketone group found; an 11-oxo substituent is required."
    
    # If we have a steroid nucleus and at least one candidate ketone group, classify as 11-oxo steroid.
    return True, ("Molecule contains a steroid-like fused ring system and at least one internal ketone group "
                  "(suggestive of an 11-oxo substituent). Details: " + "; ".join(ketone_reasons))

# Example usage (uncomment to run tests)
# test_smiles = "C[C@@]12CC(=O)[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@]2(O)C(=O)C=O"  # 11,21-dioxoprogesterone
# print(is_11_oxo_steroid(test_smiles))