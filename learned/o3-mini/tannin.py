"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: tannin
Definition: 'Any of a group of astringent polyphenolic vegetable principles or compounds,
chiefly complex glucosides of catechol and pyrogallol.'

Heuristic (improved):
  – Sanitize and add explicit hydrogens so that free hydroxyl groups are visible.
  – For every fully aromatic ring, count the number of free –OH substituents 
    (i.e. an oxygen attached to the ring carbon that itself is bonded to at least one hydrogen)
    and consider the ring as a polyphenolic unit if there are at least 2.
  – For sugar rings, require a non‐aromatic ring of size 5 or 6 that contains exactly one oxygen
    (and the rest carbons) – a common pattern in glycosides.
  – Compute the molecular weight.
  
Then classify the compound as a tannin if either:
  Option A: MW is high (>=500 Da) and there are at least 2 polyphenolic units, or
  Option B: MW is at least 300 Da, and there is at least one polyphenolic aromatic ring and one sugar ring.
  
Note: This heuristic is an attempt to balance false positives/negatives.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Uses heuristic rules based on aromatic polyphenolic rings (with free –OH groups) and sugar rings.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a tannin, False otherwise
        str: Reason for classification, including details from feature counts
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule is sanitized (this computes aromaticity etc.)
    Chem.SanitizeMol(mol)
    
    # Add explicit hydrogens so that free OH groups are visible
    mol = Chem.AddHs(mol)

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Obtain ring information (each ring as a tuple of atom indices)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    aromatic_ring_count = 0
    polyphenol_unit_count = 0  # counts aromatic rings with at least 2 hydroxyl substituents
    sugar_ring_count = 0       # counts glycosidic rings
    
    # To avoid double-counting the same -OH group from two ring atoms,
    # we keep track of oxygen indices that have been counted for a given ring.
    for ring in rings:
        # Check if the ring is fully aromatic
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
            oh_in_ring = set()
            # For each atom in the ring (we only consider carbon atoms)
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue
                # Look at neighbors that are not in the ring
                for nb in atom.GetNeighbors():
                    nb_idx = nb.GetIdx()
                    if nb_idx in ring:
                        continue
                    if nb.GetAtomicNum() == 8:
                        # Check if this oxygen has at least one hydrogen attached
                        # and ensure we have not counted it already for this ring.
                        if any(neighbor.GetAtomicNum() == 1 for neighbor in nb.GetNeighbors()):
                            oh_in_ring.add(nb_idx)
            if len(oh_in_ring) >= 2:
                polyphenol_unit_count += 1
        else:
            # Consider non-aromatic rings for potential sugar rings.
            ring_size = len(ring)
            # Limiting to rings of size 5 or 6
            if ring_size not in [5, 6]:
                continue
            # Count atoms by element:
            oxygen_count = 0
            carbon_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    oxygen_count += 1
                elif atom.GetAtomicNum() == 6:
                    carbon_count += 1
            # A typical sugar ring (pyranose or furanose) has exactly one ring oxygen.
            if oxygen_count == 1 and (carbon_count == ring_size - 1):
                sugar_ring_count += 1

    reason_details = (f"MW = {mol_wt:.1f} Da, "
                      f"{aromatic_ring_count} aromatic ring(s) (of which {polyphenol_unit_count} "
                      f"have ≥2 free OH), and "
                      f"{sugar_ring_count} sugar ring(s) detected")
    
    # Heuristic criteria:
    # Option A: High MW polyphenols: at least 500 Da with ≥2 polyphenolic (OH-rich aromatic) rings
    if mol_wt >= 500 and polyphenol_unit_count >= 2:
        return True, (f"Contains {polyphenol_unit_count} polyphenolic aromatic unit(s) and high MW "
                      f"({mol_wt:.1f} Da), consistent with tannins. Details: {reason_details}")
    
    # Option B: Glycosylated tannins: at least 300 Da with at least one polyphenolic aromatic ring and one sugar ring
    if mol_wt >= 300 and polyphenol_unit_count >= 1 and sugar_ring_count >= 1:
        return True, (f"Contains {aromatic_ring_count} aromatic ring(s) (with {polyphenol_unit_count} polyphenolic unit(s)) "
                      f"and {sugar_ring_count} sugar ring(s) at MW = {mol_wt:.1f} Da, consistent with glycosylated tannins. "
                      f"Details: {reason_details}")
    
    return False, f"Criteria not met: {reason_details}"

# Example usage when running as main:
if __name__ == "__main__":
    # Test with one example: Conoideoglucoside C
    test_smiles = "O=C(O)C1=CC(OC)=C(O[C@@H]2O[C@@H]([C@@H](OC)[C@@H]([C@H]2O)O)CO)C(=C1)CC=C(C)C"
    result, reason = is_tannin(test_smiles)
    print(result, reason)