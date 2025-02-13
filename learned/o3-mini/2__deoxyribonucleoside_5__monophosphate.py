"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: 2'-deoxyribonucleoside 5'-monophosphate
Definition: A 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.
The heuristic here:
  1. Look for a furanose sugar ring: a five-membered ring that should contain exactly one oxygen atom.
  2. Within that ring, look for a carbon that connects via an exocyclic oxygen to a phosphorus atom (the phosphate at the 5'-position).
  3. Check that there is evidence for a nucleobase (at least two nitrogen atoms located outside the sugar ring).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise.
        str: A reason explaining the basis for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flag to check if we have identified an appropriate furanose ring.
    sugar_ring_found = False
    sugar_ring_atoms = set()
    phosphate_attached = False
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # Look for a five-membered ring (furanose-like)
        if len(ring) == 5:
            # Count oxygen atoms in the ring.
            ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
            if len(ring_oxygens) != 1:
                continue  # Not a typical furanose ring.
            # Mark this ring as potential sugar ring.
            sugar_ring_found = True
            sugar_ring_atoms = set(ring)
            # Look for a carbon atom in the ring that is attached via an exocyclic oxygen to a phosphorus.
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Only consider carbon atoms from the ring as sugar carbons.
                if atom.GetAtomicNum() != 6:
                    continue
                for nbr in atom.GetNeighbors():
                    # Look for an exocyclic oxygen (not part of the ring)
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in sugar_ring_atoms:
                        # Check if this oxygen is attached to a phosphorus atom.
                        for nbr2 in nbr.GetNeighbors():
                            if nbr2.GetAtomicNum() == 15:
                                phosphate_attached = True
                                break
                        if phosphate_attached:
                            break
                if phosphate_attached:
                    break
            # If both the sugar ring and phosphate attachment found, we can stop.
            if phosphate_attached:
                break
            else:
                # If this ring qualifies as a sugar ring but without a phosphate attached, continue searching.
                sugar_ring_found = sugar_ring_found or False

    # If no appropriate five-membered sugar ring is found, classification fails.
    if not sugar_ring_found:
        return False, "No five-membered (furanose) ring with one ring oxygen (expected for deoxyribose) found"
    # If sugar ring found but no phosphate connection (via exocyclic oxygen) detected, fail.
    if not phosphate_attached:
        return False, "No phosphate group attached to the sugar ring detected (expecting a 5'-phosphate linkage)"
    
    # Check for the presence of a nucleobase.
    # A nucleobase should be attached outside the sugar ring.
    # We use a simple heuristic: at least two nitrogen atoms outside the sugar ring.
    nitrogen_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetIdx() not in sugar_ring_atoms:
            nitrogen_count += 1
    if nitrogen_count < 2:
        return False, "Insufficient evidence of a nucleobase (fewer than 2 nitrogen atoms found outside the sugar ring)"
    
    # Optionally, one could further test for aromaticity in the nucleobase region.
    # If all tests pass, classify the molecule accordingly.
    return True, "Molecule contains a furanose sugar (deoxyribose) with a 5'-phosphate attachment and evidence of a nucleobase."

# (The module ends here.)