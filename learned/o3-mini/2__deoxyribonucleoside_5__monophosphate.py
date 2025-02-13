"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: 2'-deoxyribonucleoside 5'-monophosphate
Definition: A 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.

This function uses several heuristic substructure searches:
  1. It looks for a phosphate group attached via an oxygen on a carbon (using a SMARTS pattern "COP(=O)(O)O").
  2. It searches for a five‐membered sugar ring (furanose) having exactly one ring oxygen and exactly two OH substituents on its carbon atoms.
     This pattern is consistent with 2’-deoxyribose (lacking the OH at the 2′ position).
  3. It checks that outside of that sugar-ring there is at least one aromatic heterocycle containing at least two nitrogen atoms,
     which is indicative of a nucleobase.
If all of these features are found the molecule is classified as a 2'-deoxyribonucleoside 5'-monophosphate.
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
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for a phosphate group attached via a carbon.
    # The phosphate group SMARTS here is a heuristic to look for a "COP(=O)(O)O" fragment.
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group (COP(=O)(O)O fragment) found"

    # 2. Look for a five-membered sugar ring that is consistent with deoxyribose.
    # We search for a ring of size 5 containing exactly one oxygen (the ring heteroatom).
    # Then for the carbon atoms of that ring, we count external OH groups.
    # In ribose, the sugar ring (furanose) has three hydroxyls on ring carbons while deoxyribose lacks one OH.
    ring_info = mol.GetRingInfo()
    sugar_ring_found = False
    sugar_ring_atoms = None
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            # Count the number of ring atoms that are oxygen.
            ring_oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if ring_oxygen_count != 1:
                continue  # Not a typical furanose ring.
            # For each carbon in the ring, count -OH substituents (neighbors that are oxygen with degree ==1).
            OH_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # carbon atom in the ring
                    for nbr in atom.GetNeighbors():
                        # Consider exocyclic oxygen atoms (not part of the ring) that likely represent an OH.
                        if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring and nbr.GetDegree() == 1:
                            OH_count += 1
            # For a deoxyribose ring we expect two OH groups (typically at 3' and 4', with 2' lacking OH).
            if OH_count != 2:
                continue

            # Also check that at least one carbon in this ring is attached to a phosphorus atom 
            # (indicating the phosphate is linked at the 5'-position).
            phosphate_attached = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # carbon atom
                    for nbr in atom.GetNeighbors():
                        if nbr.GetAtomicNum() == 15:  # phosphorus
                            phosphate_attached = True
                            break
                if phosphate_attached:
                    break
            
            if not phosphate_attached:
                continue
            
            # If all criteria for the sugar ring are met, we assume this is the deoxyribose part.
            sugar_ring_found = True
            sugar_ring_atoms = set(ring)
            break
    
    if not sugar_ring_found:
        return False, "No deoxyribose sugar ring with appropriate hydroxyl pattern and 5'-phosphate attachment found"
    
    # 3. Check for the presence of a nucleobase.
    # We require that there exists at least a few nitrogen atoms (at least two) that are not part of the sugar ring.
    nitrogen_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and (sugar_ring_atoms is None or atom.GetIdx() not in sugar_ring_atoms):
            nitrogen_count += 1
    if nitrogen_count < 2:
        return False, "Insufficient evidence of a nucleobase (fewer than 2 nitrogen atoms outside the sugar ring)"
    
    # Optionally, one might further check for an aromatic ring containing at least two N atoms.
    # Here we use a simple heuristic based on overall nitrogen count.
    
    # If all tests pass, we classify the molecule as a 2'-deoxyribonucleoside 5'-monophosphate.
    return True, "Molecule contains a deoxyribose sugar with a 5'-phosphate and a nucleobase."
    
# (The module ends here.)