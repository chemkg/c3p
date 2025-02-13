"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: any organic heteromonocyclic compound with a structure based on a dihydropyrrole (pyrroline)
Definition: A pyrroline derivative has a five‐membered ring containing exactly one nitrogen and four carbons,
and inside the ring there is exactly one double bond. Because many molecules are fused or contain ambiguous bond orders,
this function forces Kekulization (to reveal explicit bond orders), rejects rings with any aromatic bonds,
and (if the molecule is large and the candidate ring is fused with others) rejects the candidate.
Note that this approach is heuristic and may need further refinement.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline derivative (i.e. based on a dihydropyrrole ring)
    by checking for a five‐membered ring that contains exactly one nitrogen and four carbons,
    with exactly one double bond that is internal to the ring.
    
    Additional heuristic filters are applied: if any bond in the candidate ring was aromatic in the original
    molecule (suggesting a fully conjugated ring) the ring is rejected, and if the candidate ring is fused (i.e.
    its atoms appear in more than one ring) in a very large molecule then the ring is less likely to be the core motif.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a suitable dihydropyrrole ring is found, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Require that molecule is organic (has at least one carbon atom).
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule lacks carbon atoms; not organic."
    
    # Remember a copy of the original aromatic flags (used for filtering later)
    orig_arom = {}
    for bond in mol.GetBonds():
        orig_arom[bond.GetIdx()] = bond.GetIsAromatic()
    
    # Force Kekulization so that explicit bond orders show up.
    try:
        # clearAromaticFlags=True converts aromatic bonds into alternating single/double bonds if possible.
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        # If kekulization fails, then we can assume the aromaticity is ambiguous.
        return False, "Kekulization failed; ambiguous aromaticity."

    # Get ring information.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # For a quick size check later.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Iterate over each ring candidate.
    for ring in ring_info:
        # Consider only five-membered rings.
        if len(ring) != 5:
            continue
        
        # Count numbers of N and C atoms in the ring.
        n_count = 0
        c_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                n_count += 1
            elif atom.GetAtomicNum() == 6:
                c_count += 1
        # Look only at rings that have exactly 1 nitrogen and 4 carbons.
        if n_count != 1 or c_count != 4:
            continue
        
        # Now, examine consecutive bonds in the ring.
        double_bond_count = 0
        candidate_bonds = []
        aromatic_in_candidate = False
        # iterate in a cyclic manner over the ring indices
        for i in range(len(ring)):
            a1 = ring[i]
            a2 = ring[(i + 1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is None:
                continue
            candidate_bonds.append(bond)
            # If the original bond was flagged aromatic, we flag the candidate as aromatic.
            if orig_arom.get(bond.GetIdx(), False):
                aromatic_in_candidate = True
            # Count explicit (Kekulized) double bonds.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1

        # We require exactly one double bond inside the ring.
        if double_bond_count != 1:
            continue
        
        # Reject candidate if any bond in the ring was originally aromatic.
        if aromatic_in_candidate:
            continue

        # Check if the ring is fused.
        # Heuristically, if any atom in the ring is also part of another ring (besides this one)
        # then we consider the candidate fused.
        fused = False
        for idx in ring:
            count = 0
            for other in ring_info:
                if idx in other:
                    count += 1
            if count > 1:
                fused = True
                break
        # If the candidate is fused and the molecule is relatively large (suggesting a complex, polycyclic structure),
        # then we reject the candidate as a pyrroline core.
        if fused and (mol_wt > 600 or len(mol.GetAtoms()) > 40):
            continue

        # If we reach here, the candidate ring passes all tests.
        return True, "Found a five‐membered dihydropyrrole ring (1 nitrogen, 4 carbons, 1 double bond)."
    
    # If no candidate ring was found that satisfies our criteria.
    return False, "No valid five‐membered dihydropyrrole ring (1 nitrogen, 4 carbons, and 1 double bond) found."

# Example usage (for internal testing):
if __name__ == "__main__":
    test_smiles = [
        "O=C/1NCC(\\C1=C(\\O)/C=C/C(=C/C#C/C=C/C)/C)=O",  # Ravynic acid (should be true now)
        "S=C1NCCC1",  # Pyrrolidine-2-thione (should be false)
        "C1(N=CCC1)(C)C",  # 5,5-dimethyl-1-pyrroline (true)
        "O=C1N(C(O)(CC)C(=C1C(=O)C(CCCCCCC)C)O)C",  # Penicillenol D (true)
    ]
    for smi in test_smiles:
        result, reason = is_pyrroline(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*40}")