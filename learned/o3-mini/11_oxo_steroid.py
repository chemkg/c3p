"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: CHEBI/Custom: 11-oxo steroid
Definition: Any oxo steroid that has an oxo substituent at position 11.
Heuristic implementation:
  1. The SMILES is parsed.
  2. The molecule must contain an extended fused ring system resembling a steroid nucleus.
     In our approach we require either at least one five-membered ring and at least three six-membered rings,
     OR four six-membered rings (a common case for homologated steroids).
  3. The molecule must have at least one ring-bound ketone group. We require that
     (a) the ketone carbon (C=O) is part of a ring,
     (b) the double bond is to an oxygen, and
     (c) the carbonyl carbon is attached via its two remaining bonds only to other carbons that are themselves in rings.
  4. In addition, we require that the fused-ring system has a sufficient number of carbon atoms (≥15) so that the ketone
     is embedded in a steroid-like core.
Note: This heuristic does not perform full stereochemical or numbering assignment.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    Heuristic checks:
      - The code parses the molecule.
      - It checks for a steroid-like fused ring system. We accept two cases:
          (i) at least 1 five-membered ring and at least 3 six-membered rings, or
         (ii) 4 six-membered rings (to allow a homologated nucleus sometimes seen).
      - It searches for at least one ring-bound ketone group (C(=O)) where the carbonyl carbon,
         besides the oxygen double bond, is only attached to carbons that are in rings.
      - It also verifies that the carbon count of the entire ring system is at least 15.
      
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule is classified as 11-oxo steroid, False otherwise.
      str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry and ring perception are performed
    Chem.AssignStereochemistry(mol, cleanIt=True)
    AllChem.Compute2DCoords(mol)  # helps in some cases
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    count5 = sum(1 for ring in atom_rings if len(ring) == 5)
    count6 = sum(1 for ring in atom_rings if len(ring) == 6)
    
    # Require a fused ring system that is steroid-like.
    # We accept molecules with (a) at least 1 five-membered ring and at least 3 six-membered rings,
    # or (b) 4 six-membered rings (no five-membered rings).
    if count6 < 3:
        return False, f"Lacks enough six-membered rings (found {count6}, need at least 3)"
    if count5 < 1 and count6 != 4:
        return False, (f"Lacks a typical steroid nucleus: expected either (≥1 five-membered and ≥3 six-membered rings) "
                       f"or 4 six-membered rings, found {count5} five-membered and {count6} six-membered rings")
    
    # Next, collect all atoms that lie in any ring.
    ring_atom_indices = set()
    for ring in atom_rings:
        ring_atom_indices.update(ring)
    ring_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() in ring_atom_indices]
    if len(ring_carbons) < 15:
        return False, "Fused ring system has too few carbon atoms to be a steroid nucleus"
    
    # Now, search for a suitable ketone group: a C=O where the carbon is in a ring
    # and aside from the O attached via the double bond, its other neighbors are carbons that are themselves in rings.
    ketone_found = False
    ketone_details = []
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if not atom.IsInRing():
            continue
        
        # Iterate all bonds from this carbon
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() != 2:
                continue  # we look for double bonds only
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() != 8:
                continue  # must be an oxygen
            # We have a C=O bond. Now check the other neighbors of the carbon (aside from the O)
            other_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() != other.GetIdx()]
            # We want exactly two other neighbors for a trigonal (sp2) carbon.
            if len(other_neighbors) != 2:
                continue
            # Both of these neighbors should be carbons and be in rings (i.e. part of the fused nucleus)
            if all(nbr.GetAtomicNum() == 6 and nbr.IsInRing() for nbr in other_neighbors):
                ketone_found = True
                ketone_details.append(
                    f"Carbon atom {atom.GetIdx()} in a ring forms a ketone with oxygen atom {other.GetIdx()}"
                )
                break  # found one suitable ketone on this carbon; we can move to next atom
    
    if not ketone_found:
        return False, "No ring-bound ketone group (with both other substituents as ring carbons) found – an 11-oxo is required."
    
    # If we have a fused steroid-like ring system AND a candidate ketone, assume the molecule is an 11-oxo steroid.
    reason = ("Molecule contains a steroid-like fused ring system ("
              f"{count5} five-membered and {count6} six-membered rings, with {len(ring_carbons)} ring carbons) and at least one "
              "internal ketone group (suggestive of an 11-oxo substituent). Details: " + "; ".join(ketone_details))
    return True, reason

# Example usage (for testing):
if __name__ == "__main__":
    # a few test examples; you can add more SMILES strings from the provided list.
    test_smiles_list = [
        "O=C1C=C2C=CC3=C4[C@]5([C@@H]([C@H](C)[C@@H](C5)[C@@H](O)[C@@H](C(C)C)C)CC4)CC([C@]3([C@@]2(C)CC1)O)=O",  # Emesterone B
        "[H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(C)CC(=O)[C@@]1([H])[C@@]2([H])C[C@H](O)C2=CC(=O)C=C[C@]12C"  # 6alpha-hydroxyprednisone
    ]
    for smi in test_smiles_list:
        result, detail = is_11_oxo_steroid(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Detail:", detail)
        print("--------------------------------------------------")