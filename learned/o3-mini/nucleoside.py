"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: Nucleoside
Definition:
  “A nucleoside is an N-glycosyl compound that has both a nucleobase – normally adenine, guanine,
   xanthine, thymine, cytosine or uracil – and either a ribose or deoxyribose as functional parents.”
Improved strategy:
  1. Parse the SMILES string.
  2. Exclude molecules with phosphorus (e.g. phosphate groups => nucleotide).
  3. Identify candidate sugar rings by iterating over all rings of size 5 that contain exactly 4 carbons and 1 oxygen.
     Additionally, require that each atom in that ring belongs to only one ring (to avoid fused systems).
  4. Identify candidate nucleobase substructures using canonical SMARTS patterns, and if none are found, fall back
     to a heuristic of rings (size 5 or 6) with at least 2 nitrogen atoms.
  5. Look for a single bond from any carbon in a sugar ring to a nitrogen that is part of the nucleobase.
  
Returns:
   (True, <reason>) if the SMILES is classified as a nucleoside;
   (False, <reason>) otherwise.
   If parsing fails, returns (False, "Invalid SMILES string").
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # -------------------------------
    # Step 1. Exclude molecules with phosphorus (likely nucleotides, not nucleosides)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Phosphate group detected; molecule is likely a nucleotide, not a nucleoside"
    
    # -------------------------------
    # Step 2. Identify candidate sugar rings.
    # For a ribose or deoxyribose ring, we expect a five-membered ring with exactly 4 carbons and 1 oxygen.
    # Additionally, each atom in a genuine sugar should not be fused into another ring.
    ring_info = mol.GetRingInfo().AtomRings()
    # Count how many rings each atom belongs to.
    atom_ring_counts = {i:0 for i in range(mol.GetNumAtoms())}
    for ring in ring_info:
        for idx in ring:
            atom_ring_counts[idx] += 1

    sugar_rings = []   # will store sets of atom indices that make up candidate sugar rings
    for ring in ring_info:
        if len(ring) != 5:
            continue
        oxygen_count = 0
        carbon_count = 0
        valid_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            a_num = atom.GetAtomicNum()
            # In a ribose ring, only carbons (6) and one oxygen (8) are expected.
            if a_num == 8:
                oxygen_count += 1
            elif a_num == 6:
                carbon_count += 1
            else:
                valid_ring = False
                break
            # To avoid fused systems, require each atom in the ring belongs to only one ring.
            if atom_ring_counts[idx] > 1:
                valid_ring = False
                break
        if valid_ring and oxygen_count == 1 and carbon_count == 4:
            sugar_rings.append(set(ring))

    if not sugar_rings:
        return False, "No sugar moiety (furanose ring with 4 carbons and 1 oxygen) detected"

    # -------------------------------
    # Step 3. Identify candidate nucleobase regions.
    # First, try canonical SMARTS patterns for adenine, guanine, cytosine, thymine/uracil, and xanthine.
    nucleobase_atoms = set()
    canonical_smarts = {
        'adenine': "c1ncnc2ncnc12",
        'guanine':  "c1nc2c(n1)[nH]cnc2=O",
        'cytosine': "n1c(=O)ncn1",
        # thymine and uracil share the core pattern (thymine has one extra methyl group)
        'thymine/uracil': "n1c(=O)[nH]c(=O)n1",
        'xanthine': "c1[nH]c(=O)nc(=O)n1"
    }
    for name, smarts in canonical_smarts.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            nucleobase_atoms.update(match)

    # Fallback: If canonical patterns yield nothing, use a heuristic:
    # look for any ring (5- or 6-membered) with at least 2 nitrogen atoms.
    if not nucleobase_atoms:
        for ring in ring_info:
            if len(ring) not in (5, 6):
                continue
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_count >= 2:
                nucleobase_atoms.update(ring)

    if not nucleobase_atoms:
        return False, "No nucleobase detected via canonical SMARTS or fallback heterocycle patterns"

    # -------------------------------
    # Step 4. Identify the N-glycosidic bond connecting the sugar and the nucleobase.
    # For each candidate sugar ring, examine each carbon atom in the ring.
    # The sugar carbon should be bonded (by a single bond) to a nitrogen that is part of the nucleobase region.
    for sugar in sugar_rings:
        for idx in sugar:
            atom = mol.GetAtomWithIdx(idx)
            # We only check the carbons in the sugar ring.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # We want to leave the sugar ring.
                if nbr_idx in sugar:
                    continue
                if nbr.GetAtomicNum() != 7:  # must be a nitrogen
                    continue
                # Only consider it if the nitrogen is part of the nucleobase (i.e. in the match).
                if nbr_idx not in nucleobase_atoms:
                    continue
                bond = mol.GetBondBetweenAtoms(idx, nbr_idx)
                if bond is None:
                    continue
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                # Found a single bond connecting a carbon in a sugar ring and a nitrogen in the nucleobase.
                return True, "Molecule contains an N-glycosidic bond connecting a sugar moiety to a nucleobase"
    
    return False, "Found sugar and nucleobase motifs but no connecting N-glycosidic bond detected"

# Example usage:
if __name__ == "__main__":
    # Test with a canonical nucleoside (adenosine)
    test_smiles = "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"
    result, reason = is_nucleoside(test_smiles)
    print(result, reason)