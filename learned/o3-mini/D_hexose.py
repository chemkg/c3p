"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: D-hexose
A D-hexose is defined as a hexose (a monosaccharide containing 6 carbons) 
in which the sugar ring (either pyranose or furanose) bears an exocyclic CH2OH group 
attached to a ring carbon that has the “R” CIP configuration (i.e. D configuration at 
the stereocenter that in common nomenclature is called position 5).
Note: This heuristic only accepts molecules that appear to be a single hexose unit.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule (given as SMILES) is a D-hexose.
    
    The procedure is as follows:
    1. Parse the SMILES and add explicit hydrogens so that hydrogen counts are reliable.
    2. Check that the molecule contains exactly 6 carbon atoms (as expected for a free hexose).
    3. Find all rings of interest (size 5 or 6) that have exactly one ring oxygen (typical for sugar rings).
    4. For each such ring, examine each ring carbon:
         • Look for an exocyclic neighbor (i.e. not in the ring) that is carbon and that 
           has exactly two hydrogens (a “CH2” part) and that itself is attached to at least one oxygen 
           of degree 1 (representing the –OH group).
         • If such an exocyclic substituent is found, retrieve the CIP code of the ring carbon.
           For a D-hexose the chiral center bearing the –CH2OH (position 5) should be "R".
    5. If any candidate is found that meets the above and the molecule as a whole has exactly 6 carbons,
       return True along with an explanation.
    
    Args:
        smiles (str): input SMILES string
        
    Returns:
        bool: True if molecule is classified as a D-hexose, False otherwise
        str: Explanation of the classification result.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so CH counts are correct
    mol = Chem.AddHs(mol)
    # Assign stereochemistry (this sets _CIPCode properties)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # First, as an overall filter for free hexoses, count total carbons.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons != 6:
        return False, f"Total carbon count is {total_carbons} (expected 6 for a hexose)"
    
    # Get ring information in the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    candidate_reasons = []
    # Loop over rings of interest: pyranose (6-membered) or furanose (5-membered)
    for ring in atom_rings:
        if len(ring) not in [5, 6]:
            continue  # skip rings that are not typical sugar rings
            
        # In a valid sugar ring we expect exactly one ring oxygen.
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(ring_oxygens) != 1:
            continue
        
        # For each carbon in the ring, check if it bears an exocyclic CH2OH substituent.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # only consider carbon atoms
            # In many hexoses the carbon bearing the exocyclic CH2OH group is at the non-anomeric position.
            # Look at neighbors that are not in the ring.
            exo_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring]
            for nbr in exo_neighbors:
                if nbr.GetAtomicNum() != 6:
                    continue  # we want a carbon substituent
                # Check that this substituent carbon is CH2 (i.e. exactly 2 hydrogens attached)
                # (GetTotalNumHs already includes both explicit and implicit Hs.)
                if nbr.GetTotalNumHs() != 2:
                    continue
                # Now check that this carbon has at least one attached oxygen that is terminal (degree 1)
                oh_found = False
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetAtomicNum() == 8 and subnbr.GetDegree() == 1:
                        oh_found = True
                        break
                if not oh_found:
                    continue
                # We have found an exocyclic CH2OH group attached to the ring carbon at atom idx.
                cip = atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None
                if cip == "R":
                    return True, (f"Found sugar ring (ring of size {len(ring)}) with a CH2OH group on ring carbon "
                                  f"(atom idx {idx}) that has CIP code 'R' (D-hexose).")
                else:
                    candidate_reasons.append(
                        f"Found CH2OH on ring carbon (atom idx {idx}) but its CIP code is {cip if cip else 'undefined'}"
                    )
    # End loop over rings and atoms

    if candidate_reasons:
        return False, "Candidate sugar ring(s) found but none had the expected R configuration at the CH2OH-bearing carbon: " + "; ".join(candidate_reasons)
    
    return False, "No valid sugar ring with an exocyclic CH2OH group having the expected R configuration was found."

# Example usage for testing:
if __name__ == "__main__":
    test_examples = [
        # True positives
        ("O1[C@@H]([C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O)CO", "beta-D-altropyranose"),
        ("OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O", "beta-D-gulose"),
        ("OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O", "D-gulopyranose"),
        # False negatives (if the ring pattern isn’t recognized)
        ("O1[C@@H]([C@H](O)[C@@H](O)[C@H]1O)[C@H](O)CO", "alpha-D-glucofuranose"),
    ]
    for smi, name in test_examples:
        valid, reason = is_D_hexose(smi)
        print(f"{name}: {valid} -- {reason}")