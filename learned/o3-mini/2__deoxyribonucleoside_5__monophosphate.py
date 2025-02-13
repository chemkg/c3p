"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: 2'-deoxyribonucleoside 5'-monophosphate
Definition: A 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.

Heuristic:
  1. Identify a candidate deoxyribose sugar: a five-membered ring containing exactly one oxygen and four carbons.
  2. For each candidate ring, look for exocyclic oxygen(s) attached to its carbon atoms. For each such oxygen,
     check if it is linked to a phosphorus atom. Count the number of such sugar–O–P connections.
  3. Accept the candidate only if it shows exactly one phosphate connection (i.e. a 5'-monophosphate, not cyclic).
  4. Verify that the molecule contains at least two nitrogen atoms outside the candidate sugar ring (evidence of a nucleobase).
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
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # We will search for a candidate ring with exactly 5 atoms (4 carbons and 1 oxygen).
    candidate_ring = None
    candidate_attachment_count = 0  # Count of sugar -> exocyclic O -> P links
    candidate_ring_set = set()
    
    for ring in rings:
        if len(ring) != 5:
            continue  # Only consider 5-membered rings
        
        oxygen_count = 0
        carbon_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                carbon_count += 1
        # A deoxyribose ring should have exactly 1 oxygen and 4 carbons.
        if oxygen_count != 1 or carbon_count != 4:
            continue
        
        # Now check for phosphate attachment.
        # We expect a phosphate (phosphorus atom, Z=15) to be attached via an exocyclic oxygen.
        attachment_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We only consider carbon atoms in the ring for phosphate connectivity.
            if atom.GetAtomicNum() != 6:
                continue
            # Look at neighbors of the sugar carbon that are NOT in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # We expect an oxygen bridging the sugar carbon to a phosphate
                if nbr.GetAtomicNum() == 8:
                    # Look at neighbors of the oxygen (other than the sugar carbon)
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() == atom.GetIdx():
                            continue
                        if nbr2.GetAtomicNum() == 15:
                            attachment_count += 1
                            # We break out of the inner loop once we find a phosphate for this branch.
                            break
        # If exactly one phosphate connection is found, we choose this ring candidate.
        if attachment_count == 1:
            candidate_ring = ring
            candidate_attachment_count = attachment_count
            candidate_ring_set = set(ring)
            break
        # Otherwise, if a ring is found but with attachment_count==0 or >1, skip it.
    
    # If no candidate ring was found, then fail.
    if candidate_ring is None:
        # If we found a sugar-like ring but no proper monophosphate attachment, we give a specific message.
        if any(len(r)==5 and sum(1 for idx in r if mol.GetAtomWithIdx(idx).GetAtomicNum()==8)==1 and 
               sum(1 for idx in r if mol.GetAtomWithIdx(idx).GetAtomicNum()==6)==4 for r in rings):
            return False, "No proper 5'-phosphate attachment detected on the candidate sugar ring (either none or multiple attachments found)."
        return False, "No five-membered deoxyribose-like ring (4 carbons + 1 oxygen) detected."
    
    # Now verify nucleobase evidence:
    # Count nitrogen atoms that are not part of the candidate sugar ring.
    nucleobase_nitrogens = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetIdx() not in candidate_ring_set:
            nucleobase_nitrogens += 1
    if nucleobase_nitrogens < 2:
        return False, "Insufficient evidence of a nucleobase (fewer than 2 nitrogen atoms found outside the sugar ring)."
    
    # Optionally, check the molecular weight (nucleotides tend to be >300 Da) 
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight is too low to be a typical deoxyribonucleoside monophosphate."
    
    # All tests passed: return positive classification.
    return True, ("Molecule contains a 2'-deoxyribose sugar (5-membered ring with 4 carbons and 1 oxygen) "
                  "with a single 5'-phosphate attachment (as indicated by a sugar–O–P connectivity) and evidence of a nucleobase.")

# If the logic is not applicable (e.g. if the structure falls far outside our heuristics)
# one could return (None, None).