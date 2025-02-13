"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: Polychlorobiphenyl
Definition: A biphenyl compound containing between 2 and 10 chlorine atoms attached 
to the two benzene rings.
Improvement: Instead of taking the first biphenyl-match and counting every Cl attached anywhere,
this version iterates over all biphenyl substructure matches. For each match the code:
  1. Ensures that the two aromatic rings are “isolated” six‐membered rings joined by a single bond.
  2. Counts only chlorine atoms that are directly attached (via a single bond) to the ring carbons 
     that form the matched biphenyl scaffold.
If any one valid match has a chlorine count between 2 and 10, the molecule is classified as a PCB.
Note:
  - Some compounds (e.g. niclofolan, Formicamycin F) contain other substituents. We choose to 
    evaluate only the Cl atoms attached to the two identified benzene rings.
  - This approach may raise false positives when other aromatic biphenyl derivatives “accidentally”
    have the right chlorine count. More elaborate rules might check for extra substituents on the 
    chosen rings.
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl is defined (for our purposes) as a biphenyl compound in which one 
    identifiable biphenyl scaffold (two six-membered benzene rings connected by a single bond) 
    has between 2 and 10 chlorine atoms directly attached onto its aromatic carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a polychlorobiphenyl, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a biphenyl scaffold: two benzene rings connected by a single bond.
    # This pattern matches two 6-membered aromatic rings connected by a single bond.
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "Molecule does not contain a biphenyl scaffold (two benzene rings connected by a single bond)."
    
    # Get all matches (each match is a tuple of atom indices covering both rings)
    all_matches = mol.GetSubstructMatches(biphenyl_pattern)
    if not all_matches:
        return False, "No valid biphenyl scaffold match found."
    
    # Try each biphenyl substructure match and check chlorine count on the scaffold.
    for match in all_matches:
        # Split the 12 atoms (there should be 6 atoms for each ring) into two rings.
        # We assume the match returns 12 indices; however, because the SMARTS pattern can match in different orders,
        # we try to extract two distinct rings of size 6.
        # (One heuristic: find two sets of 6 indices that are connected.)
        match_set = set(match)
        if len(match_set) != 12:
            continue  # unexpected, skip
        
        # Get ring info from molecule
        ri = mol.GetRingInfo()
        # Identify rings (as sets of indices) that are size 6 and a subset of the current match.
        rings_in_match = []
        for ring in ri.AtomRings():
            ring_set = set(ring)
            if len(ring_set) == 6 and ring_set.issubset(match_set):
                rings_in_match.append(ring_set)
        # We expect two distinct rings. (Sometimes the same ring may be found more than once.)
        if len(rings_in_match) < 2:
            # This match does not clearly define 2 isolated benzene rings; continue with next match.
            continue
        
        # For clarity, choose the first two distinct rings
        ring1 = rings_in_match[0]
        # Look for a ring in rings_in_match that does not equal ring1.
        ring2 = None
        for r in rings_in_match[1:]:
            if r != ring1:
                ring2 = r
                break
        if ring2 is None:
            continue
        
        # Now count chlorine atoms directly attached to atoms in the two rings.
        biphenyl_atom_indices = ring1.union(ring2)
        chlorine_neighbors = set()
        for idx in biphenyl_atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider aromatic ring carbons (atomic number 6 and in aromatic ring)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                continue
            for neighbor in atom.GetNeighbors():
                # If neighbor is not part of the biphenyl scaffold and is a chlorine, count it.
                if neighbor.GetIdx() not in biphenyl_atom_indices and neighbor.GetAtomicNum() == 17:
                    # Check that the bond is a single bond (typically Cl attaches via single bond)
                    bond = mol.GetBondBetweenAtoms(idx, neighbor.GetIdx())
                    if bond is not None and bond.GetBondTypeAsDouble() == 1.0:
                        chlorine_neighbors.add(neighbor.GetIdx())
        cl_count = len(chlorine_neighbors)
        
        # Debug: if needed, one could print the chlorine count for each match.
        # Now decide if chlorine count fits within the definition (between 2 and 10).
        if cl_count < 2:
            reason = f"Found {cl_count} chlorine substituent(s) on one biphenyl scaffold; at least 2 are required."
            # Try next match if available instead of returning immediately.
            continue
        elif cl_count > 10:
            reason = f"Found {cl_count} chlorine substituent(s) on one biphenyl scaffold; no more than 10 are allowed."
            continue
        else:
            return True, f"Molecule is a polychlorobiphenyl with {cl_count} chlorine substituent(s) on its biphenyl scaffold."
    
    # If none of the matches yielded a valid chlorine count then reject the molecule.
    return False, "No biphenyl scaffold was found with a chlorine substituent count between 2 and 10."

# For example usage:
if __name__ == "__main__":
    # Example SMILES for 2,2',3,3',5,5'-hexachlorobiphenyl:
    test_smiles = "Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)c(Cl)c1Cl"
    res, msg = is_polychlorobiphenyl(test_smiles)
    print(res, msg)