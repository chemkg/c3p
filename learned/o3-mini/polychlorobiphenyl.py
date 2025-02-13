"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: Polychlorobiphenyl
Definition: A biphenyl compound containing between 2 and 10 chlorine atoms attached 
to the two benzene rings.
Improvement:
  - Instead of simply counting chlorine atoms on any biphenyl match, this version iterates 
    through all biphenyl substructure matches.
  - For each match it extracts two distinct six‐membered aromatic rings and then checks that 
    each carbon in the biphenyl scaffold does not have any extra substituents (neighbors outside 
    the scaffold) other than chlorine. It then counts the chlorine atoms attached via a single bond.
  - If any valid match yields a chlorine count between 2 and 10, the molecule is classified accordingly.
    
Note:
  - This version purposely limits the allowed substituents on the identified biphenyl scaffold 
    to only chlorine atoms so that compounds bearing other substituents (e.g. hydroxy groups) on the 
    same ring are not mis‐classified.
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl by examining the identified biphenyl scaffold.
    A valid biphenyl scaffold must consist of two isolated six-membered aromatic rings connected by a single bond.
    Additionally, none of the ring carbons may have extraneous substituents (other than chlorine).
    Then, only chlorine atoms directly attached (via single bond) to the scaffold are counted.
    If the count is between 2 and 10, the molecule is classified as a polychlorobiphenyl.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a polychlorobiphenyl, False otherwise.
        str: Explanation message.
    """
    # Parse the SMILES string into an RDKit Molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a biphenyl scaffold:
    # two benzene rings (six-membered aromatic rings) connected by a single bond.
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "Molecule does not contain a recognizable biphenyl scaffold."
    
    # Get all substructure matches for the biphenyl scaffold.
    all_matches = mol.GetSubstructMatches(biphenyl_pattern)
    if not all_matches:
        return False, "No biphenyl scaffold match found."
    
    # Iterate over all biphenyl matches.
    for match in all_matches:
        match_set = set(match)
        # Expect the match to cover 12 distinct atoms (6+6 from two benzene rings)
        if len(match_set) != 12:
            continue
        
        # Get ring information from the molecule
        ri = mol.GetRingInfo()
        rings_in_match = []
        # Identify rings within the match that are exactly size 6.
        for ring in ri.AtomRings():
            ring_set = set(ring)
            if len(ring_set) == 6 and ring_set.issubset(match_set):
                rings_in_match.append(ring_set)
        # We require two distinct six-membered rings
        if len(rings_in_match) < 2:
            continue
        ring1 = rings_in_match[0]
        ring2 = None
        for r in rings_in_match[1:]:
            if r != ring1:
                ring2 = r
                break
        if ring2 is None:
            continue
        
        # The biphenyl scaffold is the union of the two rings.
        scaffold_atoms = ring1.union(ring2)
        
        # Now check each carbon atom of the scaffold:
        # - It must be an aromatic carbon.
        # - Any neighbor outside the scaffold must be chlorine; if not, then reject this match.
        valid_scaffold = True
        chlorine_indices = set()
        for idx in scaffold_atoms:
            atom = mol.GetAtomWithIdx(idx)
            # Ensure we are only looking at aromatic carbons.
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                valid_scaffold = False
                break
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                # Only consider neighbors not in the scaffold.
                if n_idx not in scaffold_atoms:
                    # We only allow chlorine (atomic num 17) as a substituent.
                    if neighbor.GetAtomicNum() == 17:
                        # Check that the bond connecting them is a single bond.
                        bond = mol.GetBondBetweenAtoms(idx, n_idx)
                        if bond is not None and bond.GetBondTypeAsDouble() == 1.0:
                            chlorine_indices.add(n_idx)
                    else:
                        # Found an extraneous substituent: for example, OH, NO2, etc.
                        valid_scaffold = False
                        break
            if not valid_scaffold:
                break
        
        if not valid_scaffold:
            # Try next match if this scaffold has disallowed substituents.
            continue
        
        cl_count = len(chlorine_indices)
        # Check chlorine count. It must be between 2 and 10.
        if 2 <= cl_count <= 10:
            return True, f"Molecule is a polychlorobiphenyl with {cl_count} chlorine substituent(s) on its biphenyl scaffold."
        # If chlorine count is too low or too high, try next scaffold match.
    
    return False, "No biphenyl scaffold found with an allowed chlorine substituent count (2 to 10) and no extraneous substituents."

# Example usage:
if __name__ == "__main__":
    # Test with 2,2',3,3',5,5'-hexachlorobiphenyl as an example:
    test_smiles = "Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)c(Cl)c1Cl"
    result, explanation = is_polychlorobiphenyl(test_smiles)
    print(result, explanation)