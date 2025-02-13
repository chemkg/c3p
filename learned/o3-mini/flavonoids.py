"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: Flavonoids (a superclass comprising various flavonoid‐related chemotypes)
Definition (heuristic): Organic molecules whose aglycone is based on a phenyl‐substituted 
1-phenylpropane (typically a C15 or C16 skeleton) and often contains a characteristic benzopyrone 
or chalcone fragment.
Note: Because of the variability of glycosylation and other modifications the classification is heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def remove_sugars(mol):
    """
    Attempts to remove sugar rings from a molecule.
    This is done iteratively: any 5- or 6-membered ring that is non‐aromatic
    and is “isolated” (not fused to an aromatic ring) is removed.
    (This is heuristic and may over- or under-remove sugar-like fragments.)
    """
    # We loop until no more removals occur.
    mol_work = Chem.Mol(mol)
    removal_occurred = True
    while removal_occurred:
        removal_occurred = False
        ri = mol_work.GetRingInfo()
        sugar_indices = set()
        for ring in ri.AtomRings():
            ring_atoms = [mol_work.GetAtomWithIdx(idx) for idx in ring]
            n_atoms = len(ring)
            # Only consider rings that are not aromatic (sugars are aliphatic) 
            if any(atom.GetIsAromatic() for atom in ring_atoms):
                continue
            if n_atoms not in (5, 6):
                continue
            # Heuristic: if the ring atoms are mostly C and O, and at most one non C/O
            atom_syms = [atom.GetSymbol() for atom in ring_atoms]
            if all(sym in ('C','O') for sym in atom_syms):
                # Also require that none of these atoms are shared with an aromatic ring outside
                # (if fused with an aromatic ring, we expect it to be part of the flavonoid core)
                is_fused = False
                for idx in ring:
                    atom = mol_work.GetAtomWithIdx(idx)
                    for nb in atom.GetNeighbors():
                        if nb.GetIsAromatic() and nb.GetIdx() not in ring:
                            is_fused = True
                            break
                    if is_fused:
                        break
                if not is_fused:
                    sugar_indices.update(ring)
        if sugar_indices:
            # Remove atoms in descending order
            rw_mol = Chem.RWMol(mol_work)
            for idx in sorted(sugar_indices, reverse=True):
                try:
                    rw_mol.RemoveAtom(idx)
                except Exception:
                    pass
            mol_work = rw_mol.GetMol()
            Chem.SanitizeMol(mol_work, catchErrors=True)
            removal_occurred = True
    # After removal, if multiple fragments exist, keep the largest.
    frags = Chem.GetMolFrags(mol_work, asMols=True, sanitizeFrags=True)
    if frags:
        mol_work = max(frags, key=lambda m: rdMolDescriptors.CalcExactMolWt(m))
    return mol_work

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is flavonoid-like based on its SMILES string.
    
    The algorithm works as follows:
      1. Parse the SMILES string.
      2. Remove putative sugar moieties iteratively.
      3. Search for one or both of two substructure motifs:
           - A benzopyrone/flavone motif: e.g. c1ccc2oc(=O)cc2c1
           - A chalcone-like fragment: e.g. an aromatic ring connected to a carbonyl and a short chain.
      4. Compute the Murcko scaffold of the sugar-removed aglycone and check:
           - That it contains no nitrogen atoms (most flavonoids are N-free)
           - That the number of carbon atoms is in a “flavonoid‐like” range (roughly 15–22)
           - That at least two ring systems remain.
      5. If either the motif is found or the Murcko scaffold criteria are met, classify as flavonoid.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is flavonoid-like, False otherwise.
      str: Reason for classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove putative sugars from the molecule.
    aglycone = remove_sugars(mol)

    # Count aromatic rings in the aglycone.
    ri = aglycone.GetRingInfo()
    aromatic_rings = 0
    for ring in ri.AtomRings():
        if all(aglycone.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings += 1

    # Define two SMARTS patterns:
    # (a) Flavone motif (a benzopyrone): common in flavonoids.
    flavone_smarts = "c1ccc2oc(=O)cc2c1"
    flavone_pattern = Chem.MolFromSmarts(flavone_smarts)
    # (b) Chalcone-like fragment: an aromatic ring attached to a carbonyl and an sp2 carbon.
    chalcone_smarts = "c1ccc(cc1)C(=O)[C;!R]"  # [C;!R] indicates a non-ring carbon after the carbonyl.
    chalcone_pattern = Chem.MolFromSmarts(chalcone_smarts)
    
    match_flavone = aglycone.HasSubstructMatch(flavone_pattern)
    match_chalcone = aglycone.HasSubstructMatch(chalcone_pattern)
    motif_found = match_flavone or match_chalcone

    # Compute the Murcko scaffold (the “core” framework) from the aglycone.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(aglycone)
    except Exception as e:
        return False, f"Error computing Murcko scaffold: {str(e)}"

    # Count carbon and oxygen atoms in the scaffold.
    scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    scaffold_oxygens = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 8)
    # Check for nitrogen atoms – flavonoid cores are usually nitrogen-free.
    if any(atom.GetAtomicNum() == 7 for atom in scaffold.GetAtoms()):
        return False, "Aglycone scaffold contains nitrogen, which is unusual in flavonoids."

    # Count number of rings in the scaffold.
    scaffold_ri = scaffold.GetRingInfo()
    scaffold_ring_count = scaffold_ri.NumRings()

    # Use a relaxed set of conditions on the Murcko scaffold.
    scaffold_ok = (15 <= scaffold_carbons <= 22) and (scaffold_ring_count >= 2)
    
    # Formulate the reason based on our observations.
    if motif_found:
        reason = "Molecule contains a flavonoid substructure motif (flavone or chalcone-like fragment) in its aglycone."
    elif scaffold_ok:
        reason = (f"Molecule’s aglycone Murcko scaffold has {scaffold_carbons} carbons and "
                  f"{scaffold_ring_count} rings, consistent with a C15/C16 flavonoid core (allowing for modifications).")
    else:
        reason = (f"No clear flavonoid features found: motifs found = {motif_found}, "
                  f"scaffold carbons = {scaffold_carbons} (expected between 15 and 22), "
                  f"scaffold ring count = {scaffold_ring_count}.")
        return False, reason

    # Additional sanity check: require that the aglycone has at least 2 aromatic rings.
    if aromatic_rings < 2:
        return False, f"Aglycone only has {aromatic_rings} aromatic rings; flavonoids normally have at least 2."
    
    # Also check that most atoms in the aglycone are typical (C, O and small halogens)
    allowed_atomic_nums = {6, 7, 8, 9, 17, 35}  # C, N, O, F, Cl, Br (N already ruled out in scaffold)
    for atom in aglycone.GetAtoms():
        if atom.GetAtomicNum() > 20 and atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Aglycone contains atom {atom.GetSymbol()} that is unusual in flavonoids."
    
    return True, f"Flavonoid-like: {reason}"

# Example usage:
if __name__ == "__main__":
    # Test with one example (azaleatin)
    test_smiles = "COc1cc(O)cc2oc(-c3ccc(O)c(O)c3)c(O)c(=O)c12"
    result, message = is_flavonoids(test_smiles)
    print(result, message)