"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: D-hexose
A D-hexose is defined as a hexose (a monosaccharide containing 6 carbons) 
in which the sugar ring (either pyranose or furanose) bears an exocyclic CH2OH group 
attached to a ring carbon that has the “R” CIP configuration (i.e., the D configuration at 
the stereocenter that in common nomenclature is called position 5).
Note: This heuristic only accepts molecules that appear to be a single hexose unit.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose according to the heuristic:
      (1) The molecule must have exactly 6 carbons.
      (2) There must be a sugar-type ring (of size 5 or 6) containing exactly one oxygen.
      (3) One of the ring carbons must have an exocyclic substituent with a CH2OH motif.
          That is, a carbon with exactly 2 hydrogens that is bonded to a terminal OH group.
      (4) The sugar ring carbon bearing the CH2OH group must have a CIP stereochemistry of "R".
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a D-hexose, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens to ensure hydrogen counts are reliable.
    mol = Chem.AddHs(mol)
    # Force stereochemistry assignment so that _CIPCode properties are set.
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    
    # Check overall carbon count: must be 6 (free hexose).
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons != 6:
        return False, f"Total carbon count is {total_carbons} (expected 6 for a hexose)"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    candidate_reasons = []
    
    # Loop over all rings in the molecule.
    for ring in atom_rings:
        # Consider rings of size typical for sugars (5-membered furanose or 6-membered pyranose).
        if len(ring) not in [5, 6]:
            continue
            
        # In a valid sugar ring, we expect exactly one ring oxygen.
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(ring_oxygens) != 1:
            continue
            
        # Now loop through the ring atoms. We are only interested in carbons.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            
            # Look for exocyclic substituents: neighbors that are not part of the ring.
            exo_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring]
            for nbr in exo_neighbors:
                # We want the substituent to be a carbon (the CH2 part).
                if nbr.GetAtomicNum() != 6:
                    continue
                # Check if this neighboring carbon likely is CH2:
                # With explicit H’s added, a –CH2– will have exactly 2 hydrogens.
                if nbr.GetTotalNumHs() != 2:
                    continue
                    
                # Now check that the candidate exocyclic carbon is bonded to at least one terminal OH group.
                oh_found = False
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetAtomicNum() == 8 and subnbr.GetDegree() == 1:
                        oh_found = True
                        break
                if not oh_found:
                    continue
                    
                # At this point we have found an exocyclic CH2OH group attached to the ring carbon at index idx.
                # Now check for the CIP stereochemical assignment at the ring carbon.
                if not atom.HasProp('_CIPCode'):
                    candidate_reasons.append(f"Ring carbon (atom idx {idx}) lacks a CIP assignment")
                    continue
                cip = atom.GetProp('_CIPCode')
                if cip == "R":
                    return True, (f"Found sugar ring (ring of size {len(ring)}) with a CH2OH group on ring carbon "
                                  f"(atom idx {idx}) that has CIP code 'R' (D-hexose).")
                else:
                    candidate_reasons.append(
                        f"Found CH2OH on ring carbon (atom idx {idx}) but its CIP code is '{cip}'"
                    )
    
    # After examining all candidate rings.
    if candidate_reasons:
        return False, ("Candidate sugar ring(s) were found but none had the expected R configuration at the CH2OH-bearing carbon: " +
                       "; ".join(candidate_reasons))
                       
    return False, "No valid sugar ring with an exocyclic CH2OH group having the expected R configuration was found."

# Example usage for testing. Uncomment the code below to test.
if __name__ == "__main__":
    test_examples = [
        # True positives: most of the supplied examples should be recognized as D-hexoses.
        # (The list below is a subset for demonstration.)
        ("O1[C@@H]([C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O)CO", "beta-D-altropyranose"),
        ("OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O", "beta-D-gulose"),
        ("OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O", "D-gulopyranose"),
        ("O1[C@@H]([C@H](O)[C@@H](O)[C@H]1O)[C@H](O)CO", "alpha-D-glucofuranose"),
        ("[H]C(=O)[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO", "aldehydo-D-mannose"),
    ]
    for smi, name in test_examples:
        valid, reason = is_D_hexose(smi)
        print(f"{name}: {valid} -- {reason}")