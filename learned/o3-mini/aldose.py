"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: Aldose – Aldehydic parent sugars (polyhydroxy aldehydes H[CH(OH)]nC(=O)H, n>=2)
and their intramolecular hemiacetals.
Examples include beta-D-idopyranose, D-ribofuranose, L-erythrose, beta-D-glucose, etc.
"""

from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses are defined as polyhydroxy aldehydes (open-chain) or their cyclic hemiacetals.
    Typically, aldoses have between 3 and 7 carbon atoms.
    
    The algorithm handles two cases:
      1. Open-chain aldoses: The molecule is acyclic, has a terminal aldehyde group
         (SMARTS: [CX3H1](=O)) and every carbon in the backbone (traversed via DFS)
         carries at least one exocyclic hydroxyl (–OH) group.
      
      2. Cyclic (hemiacetal) aldoses: A ring of size 5 (furanose) or 6 (pyranose) is 
         identified that contains exactly one oxygen (the ring oxygen) and each ring carbon 
         bears at least one external hydroxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule appears to be an aldose, False otherwise.
        str: Detailed reason for the classification decision.
    """
    # Parse the SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms in the molecule.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    nC = len(carbons)
    if nC < 3 or nC > 7:
        return False, f"Number of carbons is {nC}, which is outside the typical range for aldoses (3-7)"
    
    # Helper function: Determines if an oxygen atom is in a hydroxyl (-OH) group,
    # i.e. it has at least one hydrogen neighbor.
    def is_hydroxyl(oxygen):
        if oxygen.GetAtomicNum() != 8:
            return False
        return any(nb.GetAtomicNum() == 1 for nb in oxygen.GetNeighbors())
    
    # Attempt 1: Open-chain aldose detection (only if molecule is acyclic).
    if mol.GetRingInfo().NumRings() == 0:
        # Use SMARTS to detect a terminal aldehyde group: a carbonyl (C=O) carbon with one hydrogen.
        aldehyde_smarts = Chem.MolFromSmarts("[CX3H1](=O)")
        aldehyde_matches = mol.GetSubstructMatches(aldehyde_smarts)
        for match in aldehyde_matches:
            aldehyde_atom = mol.GetAtomWithIdx(match[0])
            # Check that the aldehyde carbon is terminal (attached to exactly one other carbon).
            carbon_neighbors = [nb for nb in aldehyde_atom.GetNeighbors() if nb.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                continue  # Not a terminal aldehyde.
            # Traverse the carbon skeleton (backbone) via DFS.
            chain = set()
            def dfs(atom):
                if atom.GetIdx() in chain:
                    return
                chain.add(atom.GetIdx())
                for nb in atom.GetNeighbors():
                    if nb.GetAtomicNum() == 6 and nb.GetIdx() not in chain:
                        dfs(nb)
            dfs(aldehyde_atom)
            # If the DFS does not cover all carbon atoms, skip this match.
            if len(chain) != nC:
                continue
            # For every carbon in the backbone (except the aldehyde carbon),
            # check for at least one exocyclic hydroxyl group.
            open_chain_valid = True
            for idx in chain:
                if idx == aldehyde_atom.GetIdx():
                    continue
                atom = mol.GetAtomWithIdx(idx)
                found_oh = False
                # Look at neighbors not in the main carbon chain.
                for nb in atom.GetNeighbors():
                    if nb.GetIdx() in chain:
                        continue
                    if nb.GetAtomicNum() == 8 and is_hydroxyl(nb):
                        found_oh = True
                        break
                if not found_oh:
                    open_chain_valid = False
                    break
            if open_chain_valid:
                return True, ("Open-chain aldose pattern detected with terminal aldehyde "
                              "and hydroxylated carbon backbone")
    
    # Attempt 2: Cyclic hemiacetal aldose detection.
    # Look through rings of size 5 (furanose) or 6 (pyranose).
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) not in (5, 6):
            continue
        # Get the atoms in this ring.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # For a sugar ring, we expect exactly one oxygen atom (the ring oxygen).
        ring_oxygen_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
        if len(ring_oxygen_atoms) != 1:
            continue
        # The rest of the ring should be carbons.
        ring_carbon_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 6]
        if len(ring_carbon_atoms) != (len(ring) - 1):
            continue
        # For each ring carbon, check that it carries at least one exocyclic hydroxyl group.
        valid_ring = True
        for atom in ring_carbon_atoms:
            has_oh = False
            for nb in atom.GetNeighbors():
                # Only consider substituents that are not part of the ring.
                if nb.GetIdx() in ring:
                    continue
                if nb.GetAtomicNum() == 8 and is_hydroxyl(nb):
                    has_oh = True
                    break
            if not has_oh:
                valid_ring = False
                break
        if valid_ring:
            return True, ("Cyclic hemiacetal sugar pattern detected in a "
                          f"{len(ring)}-membered ring with appropriate hydroxyl decorations")
    
    # If neither open-chain nor cyclic patterns were detected, classify as non-aldose.
    return False, ("No aldose pattern detected: neither open-chain terminal aldehyde with "
                   "a fully hydroxylated backbone nor a cyclic hemiacetal sugar pattern was found")

# Example test cases (can be run if executing this file directly):
if __name__ == "__main__":
    test_cases = [
        # Open-chain aldoses:
        ("[H]C(=O)[C@@H](O)[C@H](O)CO", "aldehydo-L-ribose (open-chain)"),
        ("[H]C(=O)[C@H](O)CO", "D-glyceraldehyde (open-chain)"),
        ("[H][C@](O)(CO)[C@]([H])(O)C=O", "L-erythrose (open-chain)"),
        # Cyclic aldoses:
        ("O1[C@@H]([C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)CO", "beta-D-idopyranose (cyclic)"),
        ("O1[C@@H](O)[C@@H](O)[C@@H](O)[C@H](O)C1", "beta-D-lyxopyranose (cyclic)"),
        ("OC[C@H]1OC(O)[C@H](O)[C@@H]1O", "D-ribofuranose (cyclic)"),
        ("O1[C@H]([C@H](O)[C@H](O)C1O)[C@H](O)CO", "D-talofuranose (cyclic)"),
        ("OC[C@@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "beta-D-glucose (cyclic)")
    ]
    for smi, name in test_cases:
        result, reason = is_aldose(smi)
        print(f"SMILES: {smi} | {name} -> {result} ; Reason: {reason}")