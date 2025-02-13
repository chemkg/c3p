"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
Definition: A carbohydrate derivative in which one or more of the oxygens or hydroxy groups 
of the parent carbohydrate is replaced by sulfur or -SR, where R can be hydrogen or any group.
Heuristic:
  • For cyclic sugars: Look for a 5- or 6-membered, non-aromatic ring that would normally be a sugar (e.g. a pyranose or furanose, having one heteroatom) and that
    has at least two exocyclic hydroxyl groups. In addition, to be a thiosugar either one of the ring heteroatoms is sulfur (replacing the usual oxygen) or 
    at least one ring carbon bears an exocyclic –S substituent.
  • For acyclic sugars: Search for a sugar-like fragment which may include a carbonyl (as in ribulose) together with two or more hydroxyl groups and a nearby –S substituent.
"""

from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string using a two-route heuristic.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True and a reason if the molecule is classified as a thiosugar;
                     otherwise, False with a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so we see OH groups.
    mol = Chem.AddHs(mol)
    
    # If no sulfur at all, we cannot have a thiosugar.
    if not any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms()):
        return False, "No sulfur atom present, so not a thiosugar"
    
    # ----------------------
    # Helper functions
    # ----------------------
    def is_hydroxyl(atom, parent_idx):
        """
        Decide if an atom is an -OH substituent on a sugar ring.
        We require the atom is oxygen, carries at least one hydrogen,
        and is attached to the given parent atom.
        """
        if atom.GetSymbol() != "O":
            return False
        # Check for at least one explicit hydrogen.
        if not any(neigh.GetSymbol() == "H" for neigh in atom.GetNeighbors()):
            return False
        # Also, ensure that atom is attached to parent.
        for neigh in atom.GetNeighbors():
            if neigh.GetIdx() == parent_idx:
                return True
        return False

    def is_exocyclic_thio(atom, parent_idx):
        """
        Decide if an atom is an exocyclic -S substituent.
        For this heuristic, the substituent S must be attached to the parent
        (and not be part of a large fragment or a ring with the parent).
        """
        if atom.GetSymbol() != "S":
            return False
        # We require that the S is not itself in the same ring as parent.
        if atom.IsInRing():
            # It might be replacing the ring oxygen; we handle that separately.
            return False
        # Optionally, one can check the number of heavy neighbors:
        heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
        # Here we accept if S is bound to parent plus one small group.
        if len(heavy_neighbors) <= 2:
            return True
        return False

    # ------------------------------
    # 1. Cyclic sugar route
    # ------------------------------
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    for ring in rings:
        # Only consider typical sugar rings: 5 or 6 atoms.
        if len(ring) not in (5, 6):
            continue
        # Skip rings with any aromatic atoms.
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        # For a "sugar-like" ring, we expect that there is exactly one heteroatom
        # normally (which in a pyranose is an oxygen). Allow that this atom is S.
        hetero_atoms = [mol.GetAtomWithIdx(idx) for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
        if len(hetero_atoms) != 1:
            # If more than one heteroatom, the ring is less sugar-like.
            continue

        # Count exocyclic hydroxyl groups and exocyclic S substituents on the ring.
        oh_count = 0
        exoS_count = 0
        ring_has_S = False  # if the lone heteroatom is S rather than O.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "S":
                ring_has_S = True
            # Look at neighbors not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if is_hydroxyl(nbr, atom.GetIdx()):
                    oh_count += 1
                if is_exocyclic_thio(nbr, atom.GetIdx()):
                    exoS_count += 1

        # For a sugar ring, we require at least two -OH groups.
        if oh_count < 2:
            continue

        # To call it a thiosugar, either the ring heteroatom has been replaced with S,
        # or one or more ring carbons have an exocyclic S substituent.
        if ring_has_S or exoS_count >= 1:
            return True, ("Thiosugar identified: cyclic sugar ring (size %d) with sufficient hydroxylation " % len(ring) +
                          "and evidence of a sulfur substitution (in-ring or exocyclic).")
    # ------------------------------
    # 2. Acyclic sugar route via SMARTS
    # ------------------------------
    # Many open-chain sugars have a carbonyl (aldehyde or ketone) plus adjacent hydroxylated carbons.
    # To capture some thiosugar derivatives we require that in a sugar-like chain there is at least one
    # -S substituent near the carbonyl fragment.
    #
    # Here we use a SMARTS pattern that looks for a fragment with:
    #   • a carbon bound to S and then a carbonyl [#6]-S-[#6](=O)
    #   • followed by two chiral hydroxylated carbons ([C@H](O))
    # This is not universal but helps catch examples like 1-(methylthio)ribulose 5-phosphate.
    acyclic_smarts = Chem.MolFromSmarts("[#6]-[S]-[#6](=O)[C@H](O)[C@H](O)")
    if acyclic_smarts is not None:
        if mol.HasSubstructMatch(acyclic_smarts):
            return True, ("Thiosugar identified: acyclic sugar-like fragment containing a carbonyl, " +
                          "adjacent hydroxylated carbons and a nearby sulfur substituent." )

    # If not captured by SMARTS, try a looser search: find chains of 4-6 carbons (allowing one carbonyl) 
    # with at least two hydroxyls and at least one sulfur substituent nearby.
    atoms = mol.GetAtoms()
    visited = set()
    for atom in atoms:
        # Consider only carbon atoms not in rings.
        if atom.GetAtomicNum() != 6 or atom.IsInRing():
            continue
        # Allow both sp3 and (a single) sp2 if it is a carbonyl.
        if atom.GetHybridization().name not in ("SP3", "SP2"):
            continue
        idx = atom.GetIdx()
        if idx in visited:
            continue
        
        chain = []
        def dfs(curr, chain_set):
            chain_set.add(curr.GetIdx())
            chain.append(curr)
            for nbr in curr.GetNeighbors():
                # Restrict to carbons not in rings with SP3 character OR a single carbonyl candidate.
                if nbr.GetAtomicNum() == 6 and (not nbr.IsInRing()):
                    hyb = nbr.GetHybridization().name
                    if hyb == "SP3" or (hyb=="SP2" and any(n.GetSymbol()=="O" and n.IsDoubleBond() for n in nbr.GetNeighbors())):
                        if nbr.GetIdx() not in chain_set:
                            dfs(nbr, chain_set)
        dfs(atom, set())
        
        if not (4 <= len(chain) <= 6):
            continue
        visited.update(a.GetIdx() for a in chain)
        
        chain_oh = 0
        chain_S = 0
        # For each carbon in the chain, check substituents not in the chain.
        chain_idxs = {a.GetIdx() for a in chain}
        for ca in chain:
            for nbr in ca.GetNeighbors():
                if nbr.GetIdx() in chain_idxs:
                    continue
                if is_hydroxyl(nbr, ca.GetIdx()):
                    chain_oh += 1
                if nbr.GetSymbol() == "S":
                    chain_S += 1
        if chain_oh >= 2 and chain_S >= 1:
            return True, "Thiosugar identified: acyclic carbohydrate derivative with sufficient hydroxylation and at least one sulfur substituent nearby."
    
    return False, "No thiosugar substructure found"

# Example usage:
if __name__ == '__main__':
    # Test with one known thiosugar example.
    test_smiles = "S([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)CO)C(=NO)CC=C"
    result, reason = is_thiosugar(test_smiles)
    print(result, "->", reason)