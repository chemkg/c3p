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
  2. Identify a fused bicyclic core comprised of two 6-membered rings that share at least 2 atoms:
       • One ring must be aromatic and have at least one oxygen substituent attached via a single bond.
       • The other (saturated) ring should not be aromatic and must contain exactly one heteroatom – an oxygen.
  3. Identify a candidate attachment point on the saturated ring (a carbon that has a neighbor not in the fused core).
  4. Starting from that neighbor (which marks the beginning of the side chain), do a depth-first search 
     while not entering the fused core and while not traversing any ring atoms.
  5. In the extracted side chain, count only the carbon atoms and the carbon–carbon double bonds. 
     The chain must have between 12 and 17 carbons and either 0 or 3 double bonds.
  6. If all tests pass, classify the molecule as a tocol.
  
Note: This algorithm has been modified to allow some heteroatom substituents (eg. hydroxyl, acetyl) in the side chain,
and to relax the requirement on the substituent oxygen on the aromatic ring (which can be esterified).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines whether the molecule specified by the SMILES string is a tocol.

    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): (True, reason) if the molecule is classified as a tocol;
                     (False, reason) otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Precompute all ring atom indices in the molecule (for DFS check later)
    ring_info = mol.GetRingInfo()
    all_ring_atoms = set()
    for ring in ring_info.AtomRings():
        all_ring_atoms.update(ring)
        
    # Helpers:
    def ring_has_ox_substituent(ring):
        """Return True if any atom in the ring has a neighbor (not in the ring) 
           that is oxygen and attached via a SINGLE bond (allowing -OH or -O-acyl)."""
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in ring:
                    continue
                # only consider single-bond oxygen substituents
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                if nb.GetAtomicNum() == 8 and bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    return True
        return False

    def ring_is_aromatic(ring):
        """Return True if all atoms in the ring are set as aromatic."""
        return all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
    
    # Search for a fused bicyclic core comprised of two 6-membered rings that share at least two atoms.
    fused_core = None
    aromatic_ring = None
    saturated_ring = None
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    for i in range(len(rings)):
        ring1 = rings[i]
        if len(ring1) != 6:
            continue
        for j in range(i+1, len(rings)):
            ring2 = rings[j]
            if len(ring2) != 6:
                continue
            shared = set(ring1).intersection(set(ring2))
            if len(shared) < 2:
                continue
            # Check if one ring is aromatic and has an oxygen substituent.
            if ring_is_aromatic(ring1) and ring_has_ox_substituent(ring1):
                aromatic_ring = set(ring1)
                saturated_ring = set(ring2)
            elif ring_is_aromatic(ring2) and ring_has_ox_substituent(ring2):
                aromatic_ring = set(ring2)
                saturated_ring = set(ring1)
            else:
                continue
            # For the saturated (non-aromatic) ring, expect exactly one non-carbon atom (should be oxygen)
            non_carbons = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in saturated_ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
            if len(non_carbons) != 1 or non_carbons[0] != "O":
                continue
            # Define the fused core as the union of the two rings.
            fused_core = aromatic_ring.union(saturated_ring)
            break
        if fused_core is not None:
            break
            
    if fused_core is None:
        return False, "Chromanol core (fused aromatic and heterocyclic rings) not found"
    
    # Now, search for an attachment point on the saturated ring.
    # The candidate is a carbon in saturated ring that has at least one neighbor outside the fused core.
    candidate_attachment = None
    for idx in saturated_ring:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue
        for nb in atom.GetNeighbors():
            if nb.GetIdx() not in fused_core:
                # Accept if neighbor is carbon (even if later extra hetero atoms may be encountered).
                candidate_attachment = (idx, nb.GetIdx())
                break
        if candidate_attachment:
            break
            
    if candidate_attachment is None:
        return False, "Side chain attachment point on the chromanol core not found"
    
    # DFS to extract the side chain. We traverse from the side-chain-start (neighbor outside the core)
    # We allow any atom (even non-carbon) as long as it is not in the fused core.
    # However, if any visited atom is part of any ring, we abort (side chain must be acyclic).
    def extract_side_chain(start_idx):
        visited = set()
        stack = [start_idx]
        side_chain = set()
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            # Check: if this atom is in any ring, then the chain is cyclic.
            if curr in all_ring_atoms:
                return None, "Side chain is cyclic (contains ring atoms)"
            side_chain.add(curr)
            curr_atom = mol.GetAtomWithIdx(curr)
            for nb in curr_atom.GetNeighbors():
                if nb.GetIdx() in fused_core:
                    continue  # do not traverse into the core
                if nb.GetIdx() not in visited:
                    stack.append(nb.GetIdx())
        return side_chain, None

    side_start = candidate_attachment[1]
    side_chain, err = extract_side_chain(side_start)
    if err:
        return False, f"Error in side chain extraction: {err}"
    
    # Count the number of carbon atoms present in the extracted side chain.
    chain_carbons = sum(1 for idx in side_chain if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if chain_carbons < 12 or chain_carbons > 17:
        return False, f"Side chain has {chain_carbons} carbons; expected between 12 and 17 for three isoprenoid units"
    
    # Count double bonds in the side chain (only count bonds between carbon atoms).
    chain_double_bonds = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in side_chain and a2 in side_chain:
            atom1 = mol.GetAtomWithIdx(a1)
            atom2 = mol.GetAtomWithIdx(a2)
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    chain_double_bonds += 1
    if chain_double_bonds not in (0, 3):
        return False, f"Side chain has {chain_double_bonds} C=C double bonds; expected 0 (saturated) or 3 (triply unsaturated)"
    
    return True, "Molecule contains a chromanol core (chroman-6-ol) with an appropriate isoprenoid side chain (tocol)"

# For testing when run as a script (example with some SMILES):
if __name__ == "__main__":
    test_smiles_list = [
        # True positives (beta-tocopherol, I(3)-Tocotrienol, etc.)
        "CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@]1(C)CCc2c(C)c(O)cc(C)c2O1",  # beta-tocopherol
        "O1C(CCC=2C1=C(C(=C(O)C2)C)C)(CCC=C(CCC=C(CCC=C(C)C)C)C)C",         # I(3)-Tocotrienol
    ]
    for s in test_smiles_list:
        result, reason = is_tocol(s)
        print(f"SMILES: {s}")
        print(f"Result: {result}, Reason: {reason}\n")