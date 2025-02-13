"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: 16beta-hydroxy steroid
A 16-hydroxy steroid in which the hydroxy group at position 16 has a beta-configuration.
We use a heuristic approach:
  1. Check that the molecule has several fused rings (at least three 6-membered rings and one 5-membered ring) 
     which is typical for a steroid nucleus (cyclopentanoperhydrophenanthrene).
  2. Among the fused rings, identify a five-membered ring. In most steroids the D-ring is five-membered.
  3. Look for a hydroxy group (-OH) attached to a carbon that is part of this five-membered ring.
  4. If the candidate carbon carries an explicit stereo marker (i.e. its chirality has been defined) we assume
     that the hydroxy at that site (position 16) is in the beta configuration.
If any of these checks fail the function returns False with a reason.
Note: This is an approximate method.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    
    A steroid typically has a fused tetracyclic core (three six-membered rings and one five-membered ring).
    In a 16beta-hydroxy steroid the -OH group attached to the C atom (presumably of the D-ring, the five-membered 
    ring) is drawn with explicit stereo markers showing beta configuration.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise.
        str: Reason for classification.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that -OH groups are visible.
    mol = Chem.AddHs(mol)
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found in molecule; not a steroid nucleus"
    
    # Count rings by size and also record rings that are likely part of the steroid nucleus.
    num_5 = 0
    num_6 = 0
    rings_5 = []  # list of rings (tuple of atom indices) with 5 atoms.
    for ring in rings:
        if len(ring) == 5:
            num_5 += 1
            rings_5.append(ring)
        elif len(ring) == 6:
            num_6 += 1
            
    # A typical steroid has 3 six-membered rings and 1 five-membered ring fused together.
    if num_6 < 3 or num_5 < 1:
        return False, "Molecule does not have the typical fused ring sizes for a steroid nucleus (3 six-membered and 1 five-membered ring)"
    
    # Now, let us attempt to find a candidate hydroxy group that is attached to an atom in one of the 5-membered rings.
    candidate_found = False
    candidate_atom_idx = None
    candidate_reason = ""
    # Loop over each atom; look for carbons that belong to a 5-membered ring.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # skip non-carbons
        idx = atom.GetIdx()
        # Check if this atom is in any 5-membered ring.
        in_5_membered = any(idx in ring for ring in rings_5)
        if not in_5_membered:
            continue
        
        # Look for an -OH group attached to this carbon:
        # That is, a neighboring oxygen that is bonded to at least one hydrogen.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Check if this oxygen has at least one hydrogen.
                # We check explicit hydrogens from the molecule.
                has_H = any(nbr2.GetAtomicNum() == 1 for nbr2 in nbr.GetNeighbors())
                if has_H:
                    # We consider this -OH as the candidate.
                    candidate_found = True
                    candidate_atom_idx = idx
                    # Check if the carbon has an explicit chiral tag.
                    chiral = atom.GetChiralTag()
                    if chiral != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                        candidate_reason = ("Found fused steroid nucleus (3 six-membered and 1 five-membered rings) "
                                            "with an -OH group attached to a carbon in the five-membered ring "
                                            "that has defined stereochemistry (assumed 16beta-hydroxy).")
                    else:
                        candidate_reason = ("Found fused steroid nucleus with an -OH group on a carbon in the five-membered ring, "
                                            "but that carbon lacks explicit stereochemistry; cannot confirm beta-configuration.")
                    break
        if candidate_found:
            break

    if candidate_found:
        # We assume that if stereochemistry is present (explicit chiral tag on candidate carbon) then the -OH is beta.
        if "defined stereochemistry" in candidate_reason:
            return True, candidate_reason
        else:
            return False, candidate_reason
    else:
        return False, "Could not find an -OH group on an atom of a 5-membered ring in the presumed steroid nucleus (position 16 candidate not found)"

# For testing purposes:
if __name__ == "__main__":
    # Examples from the prompt; you can add or remove examples as needed.
    test_smiles = [
        "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O",  # 16beta-hydroxytestosterone
        "C1=C2C(CC[C@]3([C@@]4(C[C@@H]([C@@H]([C@]4(CC[C@@]32[H])C)O)O)[H])[H])=CC(=C1)O"  # 16beta-hydroxyestradiol (example)
    ]
    for s in test_smiles:
        result, reason = is_16beta_hydroxy_steroid(s)
        print("SMILES:", s)
        print("Classification:", result)
        print("Reason:", reason)
        print("-----")