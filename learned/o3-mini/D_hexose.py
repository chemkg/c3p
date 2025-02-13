"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: D-hexose
A D-hexose is defined as a hexose (a monosaccharide containing 6 carbons)
in which the sugar ring (either pyranose or furanose) bears an exocyclic CH2OH group
attached to a ring carbon that has the “R” CIP configuration (i.e., D configuration, at position 5).
This heuristic only accepts molecules that appear to be a single hexose unit.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on the following heuristic:
      (1) The molecule must have a total of 6 carbons.
      (2) It must contain a sugar ring – either a 6-membered ring (pyranose; which will have 5 carbons and 1 oxygen)
          or a 5-membered ring (furanose; which will have 4 carbons and 1 oxygen).
      (3) One of the ring carbons must bear an exocyclic substituent which is a CH2OH group.
          For our purposes, a CH2OH group is a (non‐ring) carbon with (at least) 2 hydrogens and
          which is bonded to at least one terminal OH group (an oxygen with degree 1).
      (4) The ring carbon that carries this CH2OH group must have CIP stereochemistry of "R".
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a D-hexose, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens and assign stereochemistry.
    mol = Chem.AddHs(mol)
    # Force stereochemistry assignment so that _CIPCode properties are computed.
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    
    # Check that the total carbon count is 6.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons != 6:
        return False, f"Total carbon count is {total_carbons} (expected 6 for a hexose)"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    candidate_reasons = []
    
    # Loop over all rings that might represent a sugar ring.
    for ring in atom_rings:
        # Only consider rings of the sizes typically found in sugar rings.
        if len(ring) not in (5, 6):
            continue
            
        # Within the ring, we expect exactly one heteroatom oxygen.
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(ring_oxygens) != 1:
            continue
            
        # For each atom in the ring that is a carbon:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            
            # Look for exocyclic neighbor(s); that is, neighbors not in the ring.
            exo_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring]
            for nbr in exo_neighbors:
                if nbr.GetAtomicNum() != 6:
                    continue
                # Instead of requiring exactly 2 hydrogens (which may be too strict in some cases),
                # we require at least 2 attached hydrogens.
                if nbr.GetTotalNumHs() < 2:
                    continue
                    
                # Check that this exocyclic carbon is attached to at least one terminal hydroxyl group.
                oh_found = False
                for subnbr in nbr.GetNeighbors():
                    # Look for an oxygen that is terminal (i.e., degree 1)
                    if subnbr.GetAtomicNum() == 8 and subnbr.GetDegree() == 1:
                        oh_found = True
                        break
                if not oh_found:
                    continue
                
                # We have found a candidate CH2OH group attached to the ring carbon at idx.
                # Check the stereochemical configuration at the ring carbon.
                if atom.HasProp('_CIPCode'):
                    cip = atom.GetProp('_CIPCode')
                    if cip == "R":
                        return True, (f"Found sugar ring (size {len(ring)}) with a CH2OH group on ring carbon (atom idx {idx}) "
                                      "which has CIP code 'R' (D-hexose).")
                    else:
                        candidate_reasons.append(
                            f"Ring carbon (atom idx {idx}) bearing CH2OH has CIP code '{cip}' instead of 'R'."
                        )
                else:
                    # If no CIP code available, we try to inspect the chiral tag.
                    chiral_tag = atom.GetChiralTag()
                    if chiral_tag == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                        # RDKit uses CW for R in some cases.
                        return True, (f"Found sugar ring (size {len(ring)}) with a CH2OH group on ring carbon (atom idx {idx}) "
                                      "which has chiral tag CHI_TETRAHEDRAL_CW (interpreted as 'R').")
                    else:
                        candidate_reasons.append(
                            f"Ring carbon (atom idx {idx}) bearing CH2OH lacks a recognized CIP assignment (chiral tag: {chiral_tag})."
                        )
    
    # If no candidate ring was found that satisfies the criteria.
    if candidate_reasons:
        reason = ("Candidate sugar ring(s) were found but none had the expected 'R' configuration at the CH2OH-bearing carbon: " +
                  "; ".join(candidate_reasons))
        return False, reason
    else:
        return False, "No valid sugar ring with an exocyclic CH2OH group was found."

# Example usage for testing (uncomment the lines below to test the function):
# if __name__ == "__main__":
#     test_examples = [
#         # A few examples that should be recognized as D-hexoses.
#         ("O1[C@@H]([C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O)CO", "beta-D-altropyranose"),
#         ("O1[C@@]([C@H](O)[C@@H](O)[C@@H]1O)([C@H](O)CO)[H]", "beta-D-glucofuranose"),
#         ("OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O", "beta-D-gulose"),
#         ("OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O", "D-allopyranose"),
#         ("[H]C(=O)[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO", "aldehydo-D-mannose"),
#     ]
#     for smi, name in test_examples:
#         valid, reason = is_D_hexose(smi)
#         print(f"{name}: {valid} -- {reason}")