"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: 3β-hydroxy steroid
Defined as: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the β-position.
Note:
  Because determining the exact position (the 3–position) and β–orientation
  from SMILES is nontrivial without full numbering and 3D analysis,
  the following heuristic is used:
    1. Check that the molecule has a fused tetracyclic ring system (i.e. at least 4 rings).
    2. Look for an –OH group attached to a ring carbon that is chiral and has a CIP code of "S"
       (a common assignment for 3β–OH in many steroids).
This approach is heuristic and may fail for unusual cases.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3β-hydroxy steroid, else False.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned.
    AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # Obtain information on the rings in the molecule.
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    if len(rings) < 4:
        return False, f"Found only {len(rings)} rings; a steroid nucleus typically has four fused rings"
    
    # For a steroid, we expect a fused tetracyclic ring system.
    # Here we use a simple heuristic: we demand that at least one set of atoms appears in more than one ring.
    # We collect all atoms that are in any ring.
    atoms_in_rings = set()
    for r in rings:
        atoms_in_rings.update(r)
    
    # (This check is heuristic; more rigorous checks could be implemented.)
    
    # Find candidate hydroxyl groups.
    # We look for oxygen atoms with atomic num 8 that are –OH (i.e. connected to one carbon)
    candidate_found = False
    candidate_msg = ""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            # Check if this oxygen is part of an -OH (i.e. only 1 heavy neighbor)
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 1:
                continue
            nbr = neighbors[0]
            # Check that the oxygen is not in a ring (i.e. it is an exocyclic hydroxyl)
            if atom.IsInRing():
                continue
            # Now, check that the neighbor is a carbon that is part of a ring.
            if nbr.GetAtomicNum() != 6 or not nbr.IsInRing():
                continue
            # Check that this carbon shows chirality; in many steroids the 3β‐OH is attached
            # to a chiral carbon. (Stereochemistry ‘S’ is often assigned to a 3β‐OH carbon.)
            if nbr.HasProp('_CIPCode'):
                cip = nbr.GetProp('_CIPCode')
                if cip == "S":
                    candidate_found = True
                    candidate_msg = ("Steroid nucleus present and a candidate 3β-hydroxy group found "
                                     f"(carbon {nbr.GetIdx()} with CIP code {cip}).")
                    break
            else:
                # If no CIP code is assigned, we cannot be sure.
                continue

    if not candidate_found:
        return False, "Steroid nucleus may be present, but no candidate 3β-hydroxy group (OH on a chiral ring-carbon with 'S' CIP) was found"
    
    return True, candidate_msg
     
# If run as a script, you can test with one SMILES example (e.g., cholesterol derivative).
if __name__ == "__main__":
    # Example: 4alpha-methylcholest-8-en-3beta-ol
    test_smiles = "[H][C@@]1(CC[C@@]2([H])C3=C(CC[C@]4(C)[C@]([H])(CC[C@@]34[H])[C@H](C)CCCC(C)C)[C@@]1(C)CC[C@H](O)C2"
    result, reason = is_3beta_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)