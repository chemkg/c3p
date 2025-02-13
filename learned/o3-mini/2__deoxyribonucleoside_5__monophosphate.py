"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: 2'-deoxyribonucleoside 5'-monophosphate
Definition: A 2'-deoxyribonucleoside monophosphate compound with the phosphate
group in the 5'-position.
Heuristic:
  1. Identify a candidate deoxyribose sugar: a five-membered ring with exactly one oxygen (thus 4 carbons).
  2. Within that candidate ring, check that at least one of the carbon atoms is attached (exocyclic)
     via an oxygen to a phosphorus atom (i.e. the 5'-phosphate linkage).
  3. Verify that there is evidence of a nucleobase by having at least 2 nitrogen atoms outside the sugar ring.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize flags and storage for the sugar ring candidate atoms.
    sugar_ring_found = False
    phosphate_attached = False
    sugar_ring_atoms = set()

    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    # Iterate through each ring in the molecule.
    for ring in ring_info.AtomRings():
        # We want a five-membered ring.
        if len(ring) != 5:
            continue
        
        # Count the oxygen atoms and carbon atoms inside this ring.
        oxygen_count = 0
        carbon_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                carbon_count += 1
        
        # A genuine deoxyribose ring should have exactly one oxygen and four carbons.
        if oxygen_count != 1 or carbon_count != 4:
            continue

        # This ring is a candidate sugar ring.
        sugar_ring_found = True
        current_ring = set(ring)
        # Now check for the exocyclic phosphate attachment:
        # We expect a candidate sugar ring carbon to be attached to an oxygen (exocyclic: not in the ring)
        # which in turn is attached to a phosphorus atom.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Consider only carbon atoms from the sugar ring.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                # Check for an oxygen neighbor that is not part of the ring.
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in current_ring:
                    # Look at the neighbors of this oxygen.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetAtomicNum() == 15:
                            phosphate_attached = True
                            # Save the sugar ring atom indices for later nucleobase check.
                            sugar_ring_atoms = current_ring
                            break
                    if phosphate_attached:
                        break
            if phosphate_attached:
                break
        if sugar_ring_found and phosphate_attached:
            # Found a candidate ring with the 5'-phosphate linkage.
            break
        else:
            # Reset flags if this candidate ring does not have the phosphate.
            sugar_ring_found = sugar_ring_found or False
            phosphate_attached = False

    # If no suitable candidate sugar ring was found, return failure.
    if not sugar_ring_found:
        return False, "No five-membered deoxyribose-like ring (4 carbons + 1 oxygen) detected."
    if not phosphate_attached:
        return False, "No 5'-phosphate attachment detected on the sugar ring."

    # Next, verify there is evidence of a nucleobase.
    # Our heuristic is that at least two nitrogen atoms should be present outside the sugar ring.
    nucleobase_nitrogens = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetIdx() not in sugar_ring_atoms:
            nucleobase_nitrogens += 1

    if nucleobase_nitrogens < 2:
        return False, "Insufficient evidence of a nucleobase (fewer than 2 nitrogen atoms found outside the sugar ring)."

    # If all tests pass, we classify the molecule accordingly.
    return True, "Molecule contains a deoxyribose sugar (5-membered ring with 4 carbons and 1 oxygen) with a 5'-phosphate linkage and evidence of a nucleobase."

# End of module.