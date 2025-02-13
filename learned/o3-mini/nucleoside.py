"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: Nucleoside
Definition:
  “A nucleoside is an N-glycosyl compound that has both a nucleobase – normally adenine, guanine, xanthine,
   thymine, cytosine or uracil – and either a ribose or deoxyribose as functional parents.”
Improved strategy:
  1. Parse the SMILES string.
  2. Exclude molecules with phosphorus (e.g. phosphate groups => nucleotide).
  3. Identify candidate sugar rings by iterating over rings of size 5 that contain exactly 4 carbons and 1 oxygen.
     (We relax the sp3 requirement so that dehydro sugars are also accepted.)
  4. Identify nucleobase “atoms” in two ways:
       • First, try to match canonical nucleobase SMARTS patterns for adenine, guanine, cytosine, thymine/uracil, and xanthine.
       • If no match is found, use a fallback heuristic: look for rings (size 5 or 6) that contain at least 2 nitrogen atoms.
  5. Finally, check that at least one carbon atom from a sugar ring is single-bonded to a nitrogen that is part of the nucleobase.
  
Returns:
  (True, <reason>) if the SMILES is classified as a nucleoside;
  (False, <reason>) if not.
If unable to determine, it may return (None, None).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude molecules with phosphorus (e.g. phosphate groups => nucleotide)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Phosphate group detected; molecule is likely a nucleotide, not a nucleoside"
    
    # ---------------------------------------------
    # Step 1. Identify candidate sugar rings.
    # For ribose/deoxyribose, we expect a five-membered ring with exactly 4 carbons and 1 oxygen.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []
    for ring in ring_info:
        if len(ring) != 5:
            continue
        oxygen_count = 0
        carbon_count = 0
        # We relax the sp3 check so that partly unsaturated (dehydro) sugars can be considered.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:  # oxygen
                oxygen_count += 1
            elif atom.GetAtomicNum() == 6:  # carbon
                carbon_count += 1
            else:
                # If there is any unexpected atom, skip this ring.
                carbon_count = -100
                break
        if oxygen_count == 1 and carbon_count == 4:
            sugar_rings.append(set(ring))
    
    if not sugar_rings:
        return False, "No sugar moiety (furanose ring with 4 carbons and 1 oxygen) detected"
    
    # ---------------------------------------------
    # Step 2. Identify candidate nucleobase atoms.
    # First, try SMARTS patterns for canonical nucleobases.
    nucleobase_atoms = set()
    # List of SMARTS patterns; patterns are written without stereochemistry and may not match dihydro variants.
    canonical_smarts = {
        'adenine': "c1ncnc2ncnc12",
        'guanine':  "c1nc2c(n1)[nH]cnc2=O",
        'cytosine': "n1c(=O)ncn1",
        # thimine and uracil differ by a methyl on thymine; use a pattern that matches the ring core.
        'thymine/uracil': "n1c(=O)[nH]c(=O)n1",
        # xanthine has an extra carbonyl, here we allow a loose match:
        'xanthine': "c1[nH]c(=O)nc(=O)n1"
    }
    for name, smarts in canonical_smarts.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            nucleobase_atoms.update(match)
    # Fallback: if no canonical nucleobase match, look for rings (size 5 or 6) that contain at least 2 nitrogen atoms.
    if not nucleobase_atoms:
        for ring in ring_info:
            if len(ring) not in (5, 6):
                continue
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_count >= 2:
                nucleobase_atoms.update(ring)
    
    if not nucleobase_atoms:
        return False, "No nucleobase detected via canonical SMARTS or fallback heterocycle patterns"
    
    # ---------------------------------------------
    # Step 3. Check for N-glycosidic connectivity.
    # Look for a bond connecting an atom in a sugar ring (candidate anomeric carbon)
    # directly (with a single bond) to a nitrogen that belongs to the nucleobase.
    for sugar in sugar_rings:
        for idx in sugar:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # focus on carbon atoms in the sugar ring
                continue
            # Check this sugar carbon's neighbors.
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # We want the bond to go from the sugar (inside the ring) to a nucleobase atom (which should lie outside this sugar ring).
                if nbr_idx in sugar:
                    continue
                if nbr.GetAtomicNum() != 7:  # should be a nitrogen
                    continue
                # Verify that this nitrogen is part of the nucleobase pattern.
                if nbr_idx not in nucleobase_atoms:
                    continue
                bond = mol.GetBondBetweenAtoms(idx, nbr_idx)
                if bond is None:
                    continue
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                # Found a connecting N-glycosidic bond.
                return True, "Molecule contains an N-glycosidic bond connecting a sugar moiety to a nucleobase"
    
    return False, "Found sugar and nucleobase motifs but no connecting N-glycosidic bond detected"
    
# Example usage:
if __name__ == "__main__":
    # Test with one canonical nucleoside (adenosine)
    test_smiles = "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"
    result, reason = is_nucleoside(test_smiles)
    print(result, reason)