"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: tannin
Definition: 'Any of a group of astringent polyphenolic vegetable principles or compounds,
chiefly complex glucosides of catechol and pyrogallol.'

Improved heuristic:
  – Sanitize and add explicit hydrogens so that free hydroxyl groups are visible.
  – For every fully aromatic ring, count the number of free –OH substituents 
    (i.e. an oxygen neighbor that itself is bonded to at least one hydrogen).
      • In Option A (high MW tannins) the ring counts if it has ≥2 free –OH groups.
      • In Option B (glycosylated tannins) the ring counts if it has at least 1 free –OH.
  – For sugar rings, require a non‐aromatic ring of size 5 or 6 that contains exactly one oxygen (and the rest carbons).
  – Compute the molecular weight.
  
Then classify the compound as tannin if either:
  Option A: MW >= 500 Da and there are at least 2 aromatic rings with ≥2 free –OH groups, or
  Option B: MW >= 300 Da, and there is at least one aromatic ring with ≥1 free –OH group and one sugar ring.
  
This version aims to reduce false positives while flagging molecules like Conoideoglucoside C that otherwise would be missed.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string using improved heuristic rules.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a tannin, False otherwise.
        str: Detailed reason for the classification including feature counts.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize and add explicit hydrogens so that free –OH groups are visible.
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Obtain ring information
    rings = mol.GetRingInfo().AtomRings()
    
    # Counters for aromatic rings with free –OH groups:
    polyphenol_unit_count = 0   # aromatic rings having ≥2 free -OH groups
    aromatic_with_oh_count = 0  # aromatic rings having ≥1 free -OH group
    
    # Counter for sugar rings (non-aromatic rings of size 5 or 6 with exactly one oxygen)
    sugar_ring_count = 0
    
    # Process each ring
    for ring in rings:
        # Determine if the ring is fully aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            # Fully aromatic: count free –OH groups on substituents attached to ring carbons.
            free_oh_set = set()
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # We primarily consider carbon atoms in the ring.
                if atom.GetAtomicNum() != 6:
                    continue
                for nb in atom.GetNeighbors():
                    nb_idx = nb.GetIdx()
                    # Skip if the neighbor is within the ring.
                    if nb_idx in ring:
                        continue
                    # Only consider oxygen atoms.
                    if nb.GetAtomicNum() == 8:
                        # Check if the oxygen is bonded to at least one hydrogen.
                        if any(neighbor.GetAtomicNum() == 1 for neighbor in nb.GetNeighbors()):
                            free_oh_set.add(nb_idx)
            n_free_oh = len(free_oh_set)
            # Use strict count (≥2) for Option A.
            if n_free_oh >= 2:
                polyphenol_unit_count += 1
            # Use relaxed count (≥1) for Option B.
            if n_free_oh >= 1:
                aromatic_with_oh_count += 1
        else:
            # Check for sugar ring pattern among non‐aromatic rings.
            ring_size = len(ring)
            if ring_size not in [5, 6]:
                continue
            oxygen_count = 0
            carbon_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    oxygen_count += 1
                elif atom.GetAtomicNum() == 6:
                    carbon_count += 1
            # A typical sugar ring has exactly one oxygen and the rest are carbons.
            if oxygen_count == 1 and carbon_count == ring_size - 1:
                sugar_ring_count += 1

    reason_details = (f"MW = {mol_wt:.1f} Da, "
                      f"{len(rings)} total ring(s), of which {aromatic_with_oh_count} aromatic ring(s) with ≥1 free OH "
                      f"({polyphenol_unit_count} with ≥2 free OH) and {sugar_ring_count} sugar ring(s) detected")
    
    # Heuristic criteria:
    # Option A: High MW tannins require MW >= 500 and at least 2 aromatic rings with ≥2 free OH.
    if mol_wt >= 500 and polyphenol_unit_count >= 2:
        return True, (f"Contains {polyphenol_unit_count} polyphenolic aromatic unit(s) with ≥2 free OH and high MW "
                      f"({mol_wt:.1f} Da), consistent with tannins. Details: {reason_details}")
    
    # Option B: Glycosylated tannins require MW >= 300, at least one aromatic ring with ≥1 free OH, and one sugar ring.
    if mol_wt >= 300 and aromatic_with_oh_count >= 1 and sugar_ring_count >= 1:
        return True, (f"Contains {aromatic_with_oh_count} aromatic ring(s) with free OH and {sugar_ring_count} sugar ring(s) "
                      f"at MW = {mol_wt:.1f} Da, consistent with glycosylated tannins. Details: {reason_details}")
    
    return False, f"Criteria not met: {reason_details}"

# Example usage when running as main:
if __name__ == "__main__":
    # Test with one example: Conoideoglucoside C (previously missed)
    test_smiles = "O=C(O)C1=CC(OC)=C(O[C@@H]2O[C@@H]([C@@H](OC)[C@@H]([C@H]2O)O)CO)C(=C1)CC=C(C)C"
    result, reason = is_tannin(test_smiles)
    print(result, reason)