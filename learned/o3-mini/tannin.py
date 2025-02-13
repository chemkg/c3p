"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: tannin
Definition: 'Any of a group of astringent polyphenolic vegetable principles or compounds,
chiefly complex glucosides of catechol and pyrogallol.'
Heuristic:
  – Add hydrogens so that free hydroxyl groups on aromatic rings can be counted.
  – Count aromatic rings and for each count the number of –OH (free hydroxyl) substituents.
     If an aromatic ring has two or more –OH groups it is considered a polyphenolic unit.
  – Also detect “sugar‐rings” (non‐aromatic, ring size 5–7, with oxygen fraction ≥50%) to allow
    for glycoside-type tannins.
  – Then if (A) molecular weight is high (>= 500 Da) and there are at least 2 polyphenolic units,
    or (B) the molecular weight is at least 300 Da and there is at least 1 polyphenolic aromatic ring
         and one sugar ring, we classify the compound as a tannin.
Note: This heuristic is an attempt to balance the false positives and negatives shown in previous runs.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    A tannin is a polyphenolic compound that is typically either an oligomeric/polymeric
    unit of catechol/pyrogallol rings or a complex glycoside of such phenolic units.
    
    This function uses several heuristic criteria:
      1. It counts all fully aromatic rings.
      2. For each aromatic ring it counts how many free –OH groups (hydroxyl substituents)
         are attached (by adding explicit hydrogens to the molecule first).
         A ring with at least two –OH groups is considered a polyphenolic unit.
      3. It also searches for sugar-like rings (non‐aromatic rings of size 5–7 with a high oxygen fraction)
         that may indicate a glycosidic substructure.
      4. It computes the molecular weight.
    
    Then the compound is classified as a tannin if either:
       (A) the molecular weight is high (>=500 Da) and there are at least 2 polyphenolic units, or
       (B) the molecular weight is at least 300 Da, and even if only one aromatic polyphenolic unit is present,
           there is at least one sugar ring.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a tannin, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens; needed to reliably count attached hydrogens
    mol = Chem.AddHs(mol)
    
    # Get molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # tuple of tuples of atom indices
    aromatic_ring_count = 0
    polyphenol_unit_count = 0  # count aromatic rings with ≥2 free –OH groups attached
    # For counting attached –OH groups on aromatic ring atoms:
    for ring in rings:
        # Check if the ring is fully aromatic
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
            oh_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # We only consider carbon atoms in the aromatic ring
                if atom.GetAtomicNum() != 6:
                    continue
                # Look at neighbors: count if a neighbor is oxygen and it carries at least one hydrogen
                for nb in atom.GetNeighbors():
                    if nb.GetAtomicNum() == 8:
                        # In our explicit-H molecule, hydroxyl oxygen normally has at least one hydrogen.
                        if nb.GetTotalNumHs() >= 1:
                            oh_count += 1
                            # once found for this neighbor, break so as not to count the same oxygen twice
                            break
            if oh_count >= 2:
                polyphenol_unit_count += 1

    # Next: detect sugar-like rings.
    # A typical sugar ring is non-aromatic, has 5-7 atoms, and a high oxygen fraction.
    sugar_ring_count = 0
    for ring in rings:
        # Skip if ring is aromatic
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        n_atoms = len(ring)
        if n_atoms < 5 or n_atoms > 7:
            continue
        n_oxygens = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                n_oxygens += 1
        # If at least half (>=50%) of atoms in the ring are oxygens, consider this a sugar ring.
        if n_oxygens / n_atoms >= 0.5:
            sugar_ring_count += 1

    # Compose reason string details
    reason_details = (f"MW = {mol_wt:.1f} Da, "
                      f"{aromatic_ring_count} aromatic ring(s) ({polyphenol_unit_count} with ≥2 OH), "
                      f"{sugar_ring_count} sugar ring(s) detected")
    
    # Heuristic criteria:
    # Option A: Higher-molecular weight polyphenols (oligomeric tannins)
    if mol_wt >= 500 and polyphenol_unit_count >= 2:
        return True, f"Contains {polyphenol_unit_count} polyphenolic aromatic unit(s) and high MW ({mol_wt:.1f} Da), consistent with tannins. Details: {reason_details}"
    
    # Option B: Lower-molecular weight glycosylated tannins: require at least one aromatic ring plus a sugar ring.
    if mol_wt >= 300 and aromatic_ring_count >= 1 and sugar_ring_count >= 1:
        return True, f"Contains {aromatic_ring_count} aromatic ring(s) (with {polyphenol_unit_count} polyphenolic unit(s)) and {sugar_ring_count} sugar ring(s) at MW = {mol_wt:.1f} Da, consistent with glycosylated tannins. Details: {reason_details}"
    
    # Otherwise, do not classify as tannin.
    return False, f"Criteria not met: {reason_details}"

# Example usage when running as main:
if __name__ == "__main__":
    # Test with one of the provided examples: Conoideoglucoside C (previously a false negative).
    test_smiles = "O=C(O)C1=CC(OC)=C(O[C@@H]2O[C@@H]([C@@H](OC)[C@@H]([C@H]2O)O)CO)C(=C1)CC=C(C)C"
    result, reason = is_tannin(test_smiles)
    print(result, reason)