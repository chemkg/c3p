"""
Classifies: CHEBI:39437 tocol
"""
#!/usr/bin/env python3
"""
Classifies: tocol
A tocol is defined as "A chromanol with a chroman-6-ol skeleton that is substituted at position 2 by a 
saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units".

The strategy is:
  1. Parse the SMILES.
  2. Look for a fused bicyclic system comprised of two 6-membered rings:
       - One ring must be aromatic and contain at least one –OH substituent.
       - The other ring (the “saturated” heterocycle) should contain one oxygen and the remainder carbons.
  3. Identify within the saturated ring a candidate carbon (not part of the shared edge with the aromatic ring)
       that has a neighbor outside of the fused system. That neighbor should be the start of the side chain.
  4. From that attachment, do a DFS (without re-entering the core) to extract a side chain that must be:
       - Entirely acyclic (not in any ring) and composed only of carbon atoms.
       - Have between 12 and 17 carbon atoms.
       - Contain either 0 or exactly 3 double bonds.
If these criteria are met, the molecule is classified as a tocol.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines whether the molecule specified by the SMILES string is a tocol.
    
    A tocol is defined as a chromanol (chroman-6-ol) skeleton substituted at position 2 by an
    appropriate isoprenoid side chain.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): A tuple: (True, reason) if the molecule is a tocol; (False, reason) otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, get all rings in the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    # Helper: check if a ring of given atom indices is fully aromatic.
    def ring_is_aromatic(ring):
        return all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
    
    # Helper: check if an aromatic ring has an -OH substituent.
    # For each atom in the ring, check whether it has a neighbor O that is an -OH.
    def ring_has_OH(ring):
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                # skip if neighbor is in the ring already; we want a pendent OH.
                if nb.GetIdx() in ring:
                    continue
                if nb.GetAtomicNum() == 8:
                    # Check if this oxygen has a hydrogen (explicit or implicit)
                    # RDKit gives a total H count (explicit+implicit)
                    if nb.GetTotalNumHs() >= 1:
                        return True
        return False
    
    # Next, try to find a candidate fused pair of 6-membered rings.
    # We want one ring to be aromatic with an OH and one ring to be largely saturated
    # (should contain exactly one non-carbon, which we expect to be oxygen) and not aromatic.
    fused_core = None
    aromatic_ring = None
    saturated_ring = None
    for i in range(len(rings)):
        ring1 = rings[i]
        if len(ring1) != 6:
            continue
        for j in range(i+1, len(rings)):
            ring2 = rings[j]
            if len(ring2) != 6:
                continue
            # Check if the two rings are fused: they should share at least 2 atoms.
            shared = set(ring1).intersection(set(ring2))
            if len(shared) < 2:
                continue
            # Among these two rings, try to find one that is fully aromatic and has an OH.
            if ring_is_aromatic(ring1) and ring_has_OH(ring1):
                aromatic_ring = set(ring1)
                saturated_ring = set(ring2)
            elif ring_is_aromatic(ring2) and ring_has_OH(ring2):
                aromatic_ring = set(ring2)
                saturated_ring = set(ring1)
            else:
                continue
            # For the saturated ring we expect it NOT to be fully aromatic
            # and to have one oxygen.
            non_carbons = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in saturated_ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
            if len(non_carbons) != 1 or non_carbons[0] != "O":
                # not matching expected heterocycle
                continue
            # We now define the fused core as union
            fused_core = aromatic_ring.union(saturated_ring)
            break
        if fused_core is not None:
            break

    if fused_core is None:
        return False, "Chromanol core (fused aromatic and heterocyclic rings) not found"

    # Now, look for a candidate attachment point for the isoprenoid side chain.
    # We assume the side chain is attached to the saturated ring.
    # We look for a carbon atom (in the saturated ring) that has a neighbor outside the fused core.
    candidate_attachment = None
    for idx in saturated_ring:
        atom = mol.GetAtomWithIdx(idx)
        # Only consider carbon atoms in the saturated ring as possible attachment points.
        if atom.GetAtomicNum() != 6:
            continue
        for nb in atom.GetNeighbors():
            if nb.GetIdx() not in fused_core and nb.GetAtomicNum() == 6:
                candidate_attachment = (idx, nb.GetIdx())  # tuple: (attachment atom in core, neighbor atom outside)
                break
        if candidate_attachment:
            break

    if candidate_attachment is None:
        return False, "Side chain attachment point on the chromanol core not found"
    
    # Now, starting from the neighbor (outside the core) we extract the side chain via DFS.
    # We do not traverse any atoms from the fused core.
    def extract_side_chain(start_idx):
        visited = set()
        stack = [start_idx]
        side_chain = set()
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            side_chain.add(curr)
            curr_atom = mol.GetAtomWithIdx(curr)
            # We require acyclic side chain; if any side-chain atom is in a ring, it’s an error.
            if mol.GetAtomWithIdx(curr).IsInRing():
                return None, "Side chain is cyclic (contains ring atoms)"
            for nb in curr_atom.GetNeighbors():
                # Do not traverse into the fused core.
                if nb.GetIdx() in fused_core:
                    continue
                # Only allow carbon atoms.
                if nb.GetAtomicNum() != 6:
                    return None, "Side chain contains non-carbon atoms"
                if nb.GetIdx() not in visited:
                    stack.append(nb.GetIdx())
        return side_chain, None

    side_start = candidate_attachment[1]
    side_chain, err = extract_side_chain(side_start)
    if err:
        return False, f"Error in side chain extraction: {err}"
    
    # Count carbons in the side chain.
    chain_carbons = sum(1 for idx in side_chain if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if chain_carbons < 12 or chain_carbons > 17:
        return False, f"Side chain has {chain_carbons} carbons; expected between 12 and 17 for three isoprenoid units"
    
    # Count double bonds in the side chain.
    chain_double_bonds = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in side_chain and a2 in side_chain:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                chain_double_bonds += 1
    if chain_double_bonds not in (0, 3):
        return False, f"Side chain has {chain_double_bonds} double bonds; expected 0 (saturated) or 3 (triply unsaturated)"

    return True, "Molecule contains a chromanol core (chroman-6-ol) with an appropriate isoprenoid side chain (tocol)"

# For testing when run as a script (example with beta-tocopherol):
if __name__ == "__main__":
    test_smiles_list = [
        "CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@]1(C)CCc2c(C)c(O)cc(C)c2O1",  # beta-tocopherol
        "O1C(CCC=2C1=C(C(=C(O)C2)C)C)(CCC=C(CCC=C(CCC=C(C)C)C)C)C"  # I(3)-Tocotrienol
    ]
    for s in test_smiles_list:
        result, reason = is_tocol(s)
        print(f"SMILES: {s}")
        print(f"Result: {result}, Reason: {reason}\n")