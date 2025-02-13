"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: Germacranolide – A sesquiterpene lactone based on a germacrane skeleton.
This approach checks for a 5‐membered γ‐lactone ring and for a 10‐membered carbocycle 
that is (or is fused to) the lactone ring. For a candidate 10‐membered ring we require that
at least 8 of its atoms are carbons. Also, a minimal molecular weight cut‐off is applied.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on a germacrane skeleton.
    
    The classification uses the following heuristic criteria:
      1. The molecule must contain a 5‐membered lactone ring (γ‐lactone), detected using
         the SMARTS "[CX3](=O)[OX2r5]".
      2. The molecule must also have a 10‐membered ring that is composed predominantly of carbons 
         (i.e. at least 8 out of 10 atoms are carbons).
      3. To strengthen the assignment, the lactone must be fused to or adjacent to that 10‐membered ring.
      4. The overall molecular weight should be above 200 Da.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a germacranolide, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the input SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Criterion 1: Check for a 5-membered γ-lactone using a SMARTS pattern.
    lactone_smarts = "[CX3](=O)[OX2r5]"
    lactone_pat = Chem.MolFromSmarts(lactone_smarts)
    lactone_matches = mol.GetSubstructMatches(lactone_pat)
    if not lactone_matches:
        return False, "No 5-membered γ-lactone group (ring ester) found."
    
    # Criterion 2: Identify candidate 10-membered rings in the full molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings() if ring_info is not None else []
    candidate_10_rings = []
    for ring in all_rings:
        if len(ring) == 10:
            # Count number of carbon atoms in the ring.
            carbon_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if carbon_count >= 8:
                candidate_10_rings.append(set(ring))
    
    if not candidate_10_rings:
        return False, "No suitable 10-membered carbocycle (germacrane skeleton) found."
    
    # Criterion 3: Check connectivity between the lactone and the candidate 10-membered ring.
    # For each lactone match from the SMARTS, check if its atoms (or neighbors thereof)
    # are part of one of the candidate 10-membered rings.
    fused_found = False
    for lactone in lactone_matches:
        lactone_set = set(lactone)
        for ring in candidate_10_rings:
            # If the lactone match overlaps with a candidate ring, they are fused.
            if lactone_set.intersection(ring):
                fused_found = True
                break
            # Else, check if any lactone atom is bonded to any atom in the ring.
            for atom_idx in lactone_set:
                atom = mol.GetAtomWithIdx(atom_idx)
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                if set(neighbors).intersection(ring):
                    fused_found = True
                    break
                # End inner loop for a given lactone atom.
            if fused_found:
                break
        if fused_found:
            break
    if not fused_found:
        return False, "No 10-membered germacrane skeleton fused with a 5-membered lactone found."
    
    # Criterion 4: Check the overall molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 200:
        return False, "Molecular weight too low to be a sesquiterpene lactone."
    
    # If all criteria are met, classify as a germacranolide.
    return True, ("Contains a fused 10-membered germacrane skeleton and a 5-membered γ-lactone, "
                  "with a molecular weight and carbon count consistent with sesquiterpene lactones.")

# Example usage:
if __name__ == "__main__":
    # Example with Molephantinin SMILES:
    example_smiles = r"C\C=C(/C)C(=O)O[C@H]1C\C(C)=C\C(=O)\C=C(C)/[C@H](O)[C@H]2OC(=O)C(=C)[C@H]12"
    result, reason = is_germacranolide(example_smiles)
    print(result, reason)