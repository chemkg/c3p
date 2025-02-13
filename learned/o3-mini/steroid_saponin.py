"""
Classifies: CHEBI:61655 steroid saponin
"""
#!/usr/bin/env python
"""
Classifies: Steroid saponin
Definition: Any saponin derived from a hydroxysteroid.
That is, a molecule must contain (1) a steroid (hydroxysteroid) tetracyclic core
and (2) at least one sugar moiety (a glycoside ring) attached as a substituent.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    
    Our strategy (heuristic):
      1. Parse the SMILES string.
      2. Check for a steroid (hydroxysteroid) nucleus. Here we try to match a 
         simplified SMARTS for a cyclopentanoperhydrophenanthrene core with at least one hydroxyl.
         Note: This pattern is only one example and may not capture all steroids.
      3. Look for at least one sugar ring. We use a heuristic: a ring of size 5 or 6 that 
         contains exactly one ring oxygen (as in pyranoses or furanoses) and has at least two 
         exocyclic –OH groups (oxygen atoms attached as substituents to ring carbons).
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): A tuple where the first element is True if the molecule is classified as a 
                   steroid saponin and False otherwise; the second element is a message giving 
                   the reasoning.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic 1: Look for a steroid (hydroxysteroid) nucleus.
    # This is a SMARTS pattern that aims to capture a typical cyclopentanoperhydrophenanthrene core
    # with at least one hydroxyl substituent. (Note: This pattern is not perfect.)
    steroid_core_smarts = ("[$([C@@H]1CC[C@H]2[C@@H]3CC[C@]4(C)[C@@H](O)CC[C@]4(C)"
                           "[C@H]3CC[C@H]12)]")
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid (hydroxysteroid) tetracyclic core found"
    
    # Heuristic 2: Look for at least one sugar ring.
    # We assume a sugar ring is of size 5 or 6, contains exactly one ring oxygen,
    # and has at least 2 exocyclic -OH groups (i.e. oxygen atoms attached to its ring carbons that are not part of the ring).
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    sugar_found = False
    for ring in atom_rings:
        if len(ring) not in (5,6):
            continue
        # Get ring atoms
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count the number of oxygen atoms in the ring
        n_oxy_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        if n_oxy_in_ring != 1:
            continue  # typical sugar rings (pyranose or furanose) have one ring oxygen
        
        # Now count exocyclic OH substituents on ring carbons.
        oh_count = 0
        for atom in ring_atoms:
            # We are interested in carbons in the ring; sugars usually are carbons with –OH groups as substituents.
            if atom.GetAtomicNum() != 6:
                continue
            # Check neighbors that are not in the ring.
            for nbr in atom.GetNeighbors():
                # If the neighbor is oxygen and it is not a ring member, and the bond is single
                # we assume it could be a hydroxyl group.
                if nbr.GetAtomicNum() == 8 and (nbr.GetIdx() not in ring):
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        oh_count += 1
        if oh_count >= 2:
            sugar_found = True
            break

    if not sugar_found:
        return False, "No sugar ring (glycoside) found"

    # If both a steroid core and at least one sugar ring are detected, classify as steroid saponin.
    return True, "Molecule contains a hydroxysteroid nucleus with at least one sugar ring attached"

# Example usage (these can be used for testing):
if __name__ == "__main__":
    test_smiles = [
        # Jurubine (should be steroid saponin)
        "O1[C@@]2([C@@]([C@@]3([C@]([C@]4([C@@]([C@@]5([C@@](CC4)(C[C@@H](N)CC5)[H])C)(CC3)[H])[H])(C2)[H])C)([C@@H]([C@@]1(O)CC[C@H](CO[C@@H]6O[C@@H]([C@@H](O)C(O)C6O)CO)C)C)[H])[H]",
        # 17-beta-Estradiol glucuronide (has a steroid nucleus but sugar is glucuronide;
        # depending on the ring detection it may or may not pass our heuristic)
        "O([C@@H]1[C@@]2(C(C3C(CC2)C4=C(CC3)C=C(O)C=C4)CC1)C)[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)C(O)=O"
    ]
    for smi in test_smiles:
        result, reason = is_steroid_saponin(smi)
        print("SMILES:", smi)
        print("Is steroid saponin?", result)
        print("Reason:", reason)
        print("-----")