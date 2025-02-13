"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: Flavonoids (a superclass comprising entities related to flavonoid, isoflavonoid, chalcones, etc.)
Definition (approximate): Organic molecules based on derivatives of a phenyl‐substituted 1-phenylpropane possessing a C15 or C16 skeleton
or structures condensed with C6‐C3 lignan precursors.
Note: This classification is heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def remove_sugars(mol):
    """
    A very simple sugar remover.
    This function examines all rings in the molecule and (heuristically) removes rings
    that are 5 or 6 members and contain at least 50% oxygen atoms.
    After removal, only the largest fragment is kept.
    (Note: This is a heuristic and may not remove all sugar moieties.)
    """
    ri = mol.GetRingInfo()
    sugar_atom_idxs = set()
    for ring in ri.AtomRings():
        ring_size = len(ring)
        if ring_size in (5, 6):
            n_oxygen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if n_oxygen >= (ring_size/2.0):
                # mark all atoms in this ring for removal
                sugar_atom_idxs.update(ring)
    if not sugar_atom_idxs:
        return mol  # nothing to remove

    # Remove sugar atoms from a copy
    rw_mol = Chem.RWMol(mol)
    # Remove atoms in descending order of indices (so earlier removals don't change indices)
    for idx in sorted(sugar_atom_idxs, reverse=True):
        try:
            rw_mol.RemoveAtom(idx)
        except Exception:
            pass
    newmol = rw_mol.GetMol()
    # The removal may result in several disconnected fragments; keep the largest.
    frags = Chem.GetMolFrags(newmol, asMols=True, sanitizeFrags=True)
    if frags:
        newmol = max(frags, key=lambda m: rdMolDescriptors.CalcExactMolWt(m))
    return newmol

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid (or closely related) based on its SMILES string.
    
    This revised function does the following:
      - Parses the SMILES and rejects invalid molecules.
      - Removes likely sugar moieties so that the core (aglycone) is analyzed.
      - Checks that the aglycone has at least 2 aromatic rings.
      - Searches for substructure patterns typical of flavonoids (benzopyran or chalcone-like fragments).
      - Computes the Murcko scaffold of the aglycone and counts carbon atoms;
        additionally, it rejects molecules whose aglycone scaffold contains nitrogen.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is flavonoid-like, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, compute a version of the molecule with putative sugars removed.
    aglycone = remove_sugars(mol)
    
    # Check that there are at least 2 aromatic rings in the aglycone.
    ri = aglycone.GetRingInfo()
    aromatic_rings = 0
    for ring in ri.AtomRings():
        if all(aglycone.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings += 1
    if aromatic_rings < 2:
        return False, (f"Aglycone has only {aromatic_rings} aromatic ring(s); flavonoids normally have at least 2 aromatic rings.")

    # Define SMARTS patterns for a benzopyran-like motif and a chalcone-like fragment.
    # The benzopyran pattern looks for a benzene ring fused via an oxygen to a heterocycle.
    benzopyran_smarts = "c1ccc2occc2c1"
    benzopyran_pattern = Chem.MolFromSmarts(benzopyran_smarts)

    # A chalcone is roughly two aromatic rings connected by a carbonyl and a short chain.
    # This is an approximate pattern.
    chalcone_smarts = "c1ccc(cc1)C(=O)C"
    chalcone_pattern = Chem.MolFromSmarts(chalcone_smarts)

    # Check for the presence of either motif in the aglycone.
    has_flavonoid_motif = aglycone.HasSubstructMatch(benzopyran_pattern) or aglycone.HasSubstructMatch(chalcone_pattern)

    # Compute the Murcko scaffold of the aglycone.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(aglycone)
    except Exception as e:
        return False, f"Error computing Murcko scaffold: {str(e)}"

    # Count the number of carbon atoms in the scaffold.
    scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check that the scaffold does not contain nitrogen (flavonoids are typically N-free)
    if any(atom.GetAtomicNum() == 7 for atom in scaffold.GetAtoms()):
        return False, "Aglycone scaffold contains nitrogen, which is unusual in flavonoids."
    
    # Heuristic: A flavonoid is indicated if either
    #   (a) a typical flavonoid motif is found, or
    #   (b) the Murcko scaffold (from the sugar-removed molecule) has about 15 or 16 carbons.
    # (Note: The sugar removal should help lower the scaffold carbon count for glycosylated flavonoids.)
    has_correct_scaffold = scaffold_carbons in (15, 16)

    if has_flavonoid_motif:
        reason = "Molecule contains a flavonoid motif (benzopyran or chalcone-like fragment) in its aglycone."
    elif has_correct_scaffold:
        reason = f"Molecule has a Murcko scaffold (aglycone) with {scaffold_carbons} carbons, consistent with a C15/C16 flavonoid core."
    else:
        reason = ("No clear flavonoid features found: neither a benzopyran/chalcone motif was detected "
                  f"nor was a core with 15 or 16 carbons identified (scaffold carbons = {scaffold_carbons}).")
        return False, reason

    # Additional check: ensure that the molecule is largely composed of typical flavonoid atoms.
    # Flavonoids usually contain C, O (possibly small amounts of halogens); nitrogen is uncommon.
    allowed_atomic_nums = {6, 7, 8, 9, 17, 35}  # C, N, O, F, Cl, Br
    for atom in aglycone.GetAtoms():
        # If any heavy atom beyond allowed ones is present, reject.
        if atom.GetAtomicNum() > 20 and atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Aglycone contains atom {atom.GetSymbol()} that is unusual in flavonoids"
    
    return True, f"Flavonoid-like: {reason}"

# Example usage:
if __name__ == "__main__":
    # Test a few SMILES (one example is azaleatin; many more may be tested)
    test_smiles = "COc1cc(O)cc2oc(-c3ccc(O)c(O)c3)c(O)c(=O)c12"
    result, message = is_flavonoids(test_smiles)
    print(result, message)