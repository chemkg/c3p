"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: Nucleotide
A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a 
nucleoside with phosphoric acid.
This implementation attempts the following:
  1. The molecule must contain at least one phosphorus atom.
  2. Find a five‐membered (furanose‐like) ring composed of 4 carbons and 1 oxygen (the sugar).
  3. Extend the sugar set to include directly attached oxygens (e.g. the CH2OH group).
  4. Find a candidate glycosidic bond: some sugar atom bonded to a nitrogen that belongs to an aromatic ring 
     containing at least two nitrogen atoms (the nucleobase).
  5. Find at least one phosphate ester connected to a sugar–OH (by checking for an oxygen on a sugar carbon 
     that is single–bonded to a phosphorus atom, and not “acylated”).
  6. Finally, remove the atoms that belong to the nucleotide core (sugar, base, phosphate group(s)) and 
     reject if many atoms remain (i.e. the molecule contains extra “appendages” that are more typical of CoA‐derivatives).
If any step fails, a reason is returned.
Note: This approach is heuristic. There are many boundary cases and “fuzzy” structures in nucleotide chemistry.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a nucleotide, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Must contain at least one phosphorus atom
    phosphate_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphate_atoms:
        return False, "No phosphorus atom(s) found"
    
    # 2. Identify the sugar moiety: search for a 5-membered ring with exactly 1 oxygen and 4 carbons.
    ring_info = mol.GetRingInfo()
    sugar_ring = None
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            nO = sum(1 for a in atoms if a.GetSymbol() == "O")
            nC = sum(1 for a in atoms if a.GetSymbol() == "C")
            if nO == 1 and nC == 4:
                sugar_ring = set(ring)
                break
    if sugar_ring is None:
        return False, "Sugar moiety (five‐membered ring with 4 C and 1 O) not found"
    
    # 3. Extend sugar set to include exocyclic oxygens (e.g. CH2OH at the 5' position)
    sugar_extended = set(sugar_ring)
    for idx in list(sugar_ring):
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # add directly attached oxygens not already in the ring
            if nbr.GetSymbol() == "O":
                sugar_extended.add(nbr.GetIdx())
    
    # 4. Identify nucleobase: look for a glycosidic bond from a sugar atom to a nitrogen that is part of an aromatic ring.
    # We require that the ring (or one of the rings in a fused system) contains at least 2 nitrogen atoms.
    base_set = None
    for s_idx in sugar_extended:
        sugar_atom = mol.GetAtomWithIdx(s_idx)
        for nbr in sugar_atom.GetNeighbors():
            if nbr.GetIdx() in sugar_extended:
                continue
            if nbr.GetAtomicNum() != 7:
                continue
            # Check that the bond is a single bond (glycosidic link is typically single)
            bond = mol.GetBondBetweenAtoms(s_idx, nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            # Look for any ring containing this neighbor that is aromatic
            for ring in ring_info.AtomRings():
                if nbr.GetIdx() in ring:
                    ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                    # Accept even if one atom might escape aromaticity (fused systems)
                    n_arom = sum(1 for a in ring_atoms if a.GetIsAromatic())
                    if n_arom < len(ring_atoms) - 1:
                        continue
                    nN = sum(1 for a in ring_atoms if a.GetAtomicNum() == 7)
                    if nN >= 2:
                        base_set = set(ring)
                        break
            if base_set is not None:
                break
        if base_set is not None:
            break
    if base_set is None:
        return False, "Nucleobase moiety not found; no aromatic ring with ≥2 nitrogen atoms connected to sugar"
    
    # 5. Look for a phosphate ester attached to the sugar.
    # Instead of scanning all P atoms, we check for a sugar carbon (or exocyclic oxygen on sugar) that bonds to an O
    # which in turn is bonded to a phosphorus.
    phosphate_group = set()
    phosphate_found = False
    for idx in sugar_extended:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() not in ["C", "O"]:
            continue
        for nbr in atom.GetNeighbors():
            # Look at oxygens not already in sugar_ext (the hydroxyl groups on sugar)
            if nbr.GetSymbol() != "O" or nbr.GetIdx() in sugar_extended:
                continue
            # Check that this oxygen is not “acylated” (i.e. further bound to a carbonyl carbon)
            is_acyl = False
            for subnbr in nbr.GetNeighbors():
                if subnbr.GetIdx() == atom.GetIdx():
                    continue
                if subnbr.GetSymbol() == "C":
                    for b in subnbr.GetBonds():
                        # if subnbr has a double bond to an oxygen, mark as acyl
                        other = b.GetOtherAtom(subnbr)
                        if other.GetSymbol() == "O" and b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            is_acyl = True
                            break
                    if is_acyl:
                        break
            if is_acyl:
                continue
            # Now check if this oxygen is bonded to a phosphorus.
            for second_nbr in nbr.GetNeighbors():
                if second_nbr.GetSymbol() == "P":
                    # Found a phosphate ester directly attached to the sugar.
                    phosphate_found = True
                    phosphate_group.add(second_nbr.GetIdx())
                    phosphate_group.add(nbr.GetIdx())
                    # Also add other oxygens linked to the P that are part of the phosphate chain.
                    for p_nbr in second_nbr.GetNeighbors():
                        if p_nbr.GetSymbol() == "O":
                            phosphate_group.add(p_nbr.GetIdx())
                    break
            if phosphate_found:
                break
        if phosphate_found:
            break
    if not phosphate_found:
        return False, "No phosphate group found attached to a sugar hydroxyl (no condensation at 3' or 5' OH)"
    
    # 6. Build the nucleotide core from the three pieces.
    nucleotide_core = sugar_extended.union(base_set).union(phosphate_group)
    
    # Optional sanity check:
    # In a bona fide nucleotide the remainder of the molecule should be relatively small.
    all_heavy = set(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() not in [1])
    extra = all_heavy - nucleotide_core
    # Allow a few extra atoms (e.g. counterions or minor substituents) but reject if many atoms remain.
    if len(extra) > 10:
        return False, f"Extra {len(extra)} heavy atoms detected outside the nucleoside phosphate core (likely a conjugate)"
    
    # (Also one might check that the molecular weight isn’t absurdly high.)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 800:
        return False, "Molecular weight too high for a typical nucleotide"
    
    return True, "Contains nucleoside phosphate: sugar, nucleobase and phosphate connected correctly"

# (Example usage - can be commented out in production)
if __name__ == "__main__":
    # Example: propanoyl-AMP (expected to be classified as a nucleotide)
    smiles_example = "CCC(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12"
    result, reason = is_nucleotide(smiles_example)
    print(result, reason)