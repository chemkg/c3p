"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: Flavonoids (a superclass comprising various flavonoid‐related chemotypes)
Definition (heuristic): Organic molecules whose aglycone is based on derivatives of a phenyl‐substituted 
1-phenylpropane (typically a C15 or C16 skeleton, sometimes extended by condensation with C6–C3 lignans). 
Because of glycosylation and other modifications the classification is necessarily heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def remove_sugars(mol):
    """
    Attempts to remove sugar rings from a molecule.
    This is done iteratively. Any isolated (non-fused to an aromatic ring) 5- or 6-membered aliphatic
    ring (with atoms mostly C and O) is removed.
    This heuristic procedure may over– or under–remove sugar-like fragments.
    """
    mol_work = Chem.Mol(mol)
    removal_occurred = True
    while removal_occurred:
        removal_occurred = False
        ri = mol_work.GetRingInfo()
        sugar_indices = set()
        for ring in ri.AtomRings():
            ring_atoms = [mol_work.GetAtomWithIdx(idx) for idx in ring]
            n_atoms = len(ring)
            # Consider only non‐aromatic rings typically corresponding to sugars.
            if any(atom.GetIsAromatic() for atom in ring_atoms):
                continue
            if n_atoms not in (5, 6):
                continue
            # Heuristic: if all atoms are C or O, consider this as a candidate sugar ring.
            if all(atom.GetSymbol() in ('C', 'O') for atom in ring_atoms):
                # Also ensure that none is fused to an external aromatic ring.
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
            rw_mol = Chem.RWMol(mol_work)
            for idx in sorted(sugar_indices, reverse=True):
                try:
                    rw_mol.RemoveAtom(idx)
                except Exception:
                    pass
            mol_work = rw_mol.GetMol()
            # Sanitize and update before the next iteration.
            Chem.SanitizeMol(mol_work, catchErrors=True)
            removal_occurred = True
    # After removal, if multiple fragments remain, pick the largest by molecular weight.
    frags = Chem.GetMolFrags(mol_work, asMols=True, sanitizeFrags=True)
    if frags:
        mol_work = max(frags, key=lambda m: rdMolDescriptors.CalcExactMolWt(m))
    return mol_work

def is_flavonoids(smiles: str):
    """
    Determines if the molecule (given as a SMILES string) is flavonoid-like.

    The algorithm is as follows:
      1. Parse the SMILES string.
      2. Remove putative sugar moieties (common in glycosylated flavonoids).
      3. Search for typical flavonoid substructures:
             - A flavone/benzopyrone motif (e.g., "c1ccc2oc(=O)cc2c1").
             - A chalcone-like motif (an aromatic ring connected to a carbonyl and an sp2 carbon).
      4. Compute the Murcko scaffold of the sugar-removed aglycone and check:
             - That it contains no nitrogen atoms.
             - That the number of carbon atoms is roughly 15–30.
             - That it has at least two rings.
      5. Also count aromatic rings in the aglycone; many flavonoids contain at least three.
      6. The final decision is positive if either a typical substructure is found or
         if the scaffold fits the criteria and the molecule has at least three aromatic rings.
    
    Args:
      smiles (str): SMILES string for the molecule.
    
    Returns:
      (bool, str): Tuple where the first element is True if classified as flavonoid-like,
                   otherwise False; the second element gives the reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove putative sugars to expose the aglycone.
    aglycone = remove_sugars(mol)

    # Count aromatic rings in the aglycone.
    ri = aglycone.GetRingInfo()
    aromatic_rings = 0
    for ring in ri.AtomRings():
        if all(aglycone.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings += 1

    # Define SMARTS patterns for common flavonoid fragments.
    # (a) A flavone/benzopyrone motif.
    flavone_smarts = "c1ccc2oc(=O)cc2c1"
    flavone_pattern = Chem.MolFromSmarts(flavone_smarts)
    # (b) A chalcone-like motif: an aromatic ring attached to a carbonyl and then an sp2 carbon.
    chalcone_smarts = "c1ccc(cc1)C(=O)[C;!R]"  # [C;!R] indicates a non-ring carbon.
    chalcone_pattern = Chem.MolFromSmarts(chalcone_smarts)
    
    match_flavone = aglycone.HasSubstructMatch(flavone_pattern)
    match_chalcone = aglycone.HasSubstructMatch(chalcone_pattern)
    motif_found = match_flavone or match_chalcone

    # Compute the Murcko scaffold of the aglycone.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(aglycone)
    except Exception as e:
        return False, f"Error computing Murcko scaffold: {str(e)}"
    
    scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    scaffold_ring_count = scaffold.GetRingInfo().NumRings()
    
    # Exclude if any nitrogen is present in the core.
    if any(atom.GetAtomicNum() == 7 for atom in scaffold.GetAtoms()):
        return False, "Aglycone scaffold contains nitrogen, which is unusual in flavonoids."
    
    # Check for allowed atoms in the aglycone.
    allowed_atomic_nums = {6, 8, 9, 17, 35}  # C, O, F, Cl, Br (nitrogen already forbidden)
    for atom in aglycone.GetAtoms():
        an = atom.GetAtomicNum()
        if an > 20 and an not in allowed_atomic_nums:
            return False, f"Aglycone contains atom {atom.GetSymbol()} that is unusual in flavonoids."

    # Relaxed scaffold criteria: accept scaffolds with roughly 15–30 carbons and two or more rings.
    scaffold_ok = (15 <= scaffold_carbons <= 30) and (scaffold_ring_count >= 2)
    
    # Formulate our decision.
    if motif_found:
        reason = "Molecule contains a flavonoid substructure motif (flavone or chalcone-like fragment) in its aglycone."
        decision = True
    elif scaffold_ok and (aromatic_rings >= 3):
        reason = (f"Molecule’s aglycone Murcko scaffold has {scaffold_carbons} carbons and "
                  f"{scaffold_ring_count} rings with {aromatic_rings} aromatic rings, "
                  "consistent with a typical flavonoid core (allowing for modifications).")
        decision = True
    else:
        reason = (f"No clear flavonoid features found: motif_found = {motif_found}, "
                  f"scaffold carbons = {scaffold_carbons} (expected 15–30), "
                  f"scaffold rings = {scaffold_ring_count}, aromatic rings = {aromatic_rings}.")
        return False, reason

    # Additional sanity: Most flavonoid aglycones have at least three aromatic rings.
    if aromatic_rings < 3:
        return False, f"Aglycone only has {aromatic_rings} aromatic rings; typical flavonoids usually exhibit ≥3."

    return decision, f"Flavonoid-like: {reason}"

# Example usage for testing purposes:
if __name__ == "__main__":
    # Test with azaleatin – one of the example flavonoids.
    test_smiles = "COc1cc(O)cc2oc(-c3ccc(O)c(O)c3)c(O)c(=O)c12"
    result, msg = is_flavonoids(test_smiles)
    print(result, msg)