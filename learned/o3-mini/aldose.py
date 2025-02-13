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
    Aldoses are defined as polyhydroxy aldehydes (open-chain) or their cyclic hemiacetal forms.
    Typically, aldoses have between 3 and 7 carbon atoms.

    The algorithm first checks for an open-chain form:
      - It finds an aldehyde group (using [CX3H1](=O)) and ensures that the aldehyde 
        carbon is terminal (i.e. attached to exactly one other carbon).
      - The algorithm then traverses the connected carbon backbone. If all carbon atoms
        (except the aldehyde carbon) carry an exocyclic hydroxyl (–OH) group, we classify
        the molecule as an open-chain aldose.

    If no valid open-chain pattern is found, the code then looks for cyclic forms:
      - It iterates over rings of size 5 or 6.
      - For candidate rings, it requires that the ring contains exactly one oxygen (the
        ring oxygen) and that every ring carbon (i.e. the remaining atoms) has an exocyclic
        hydroxyl (–OH) group.
        
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule appears to be an aldose, False otherwise
        str: Detailed reason for the classification decision
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    nC = len(carbons)
    if nC < 3 or nC > 7:
        return False, f"Number of carbons is {nC}, which is outside the typical range for aldoses (3-7)"
    
    # Function to check if a given oxygen atom is part of a hydroxyl (-OH)
    def is_hydroxyl(oxygen):
        if oxygen.GetAtomicNum() != 8:
            return False
        # Check that it has at least one hydrogen neighbor
        h_count = sum(1 for nb in oxygen.GetNeighbors() if nb.GetAtomicNum() == 1)
        return h_count >= 1

    # Attempt 1: Open-chain aldose detection
    # Use SMARTS for an aldehyde group: a carbonyl carbon with one hydrogen.
    aldehyde_smarts = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_smarts)
    for match in aldehyde_matches:
        aldehyde_atom = mol.GetAtomWithIdx(match[0])
        # Check: aldehyde carbon should have exactly one neighbor that is carbon.
        carbon_neighs = [nb for nb in aldehyde_atom.GetNeighbors() if nb.GetAtomicNum() == 6]
        if len(carbon_neighs) != 1:
            continue  # Not a terminal aldehyde
        # Traverse the carbon backbone from the aldehyde atom using BFS.
        visited = set()
        queue = [aldehyde_atom]
        chain_atoms = []
        while queue:
            atom = queue.pop()
            if atom.GetAtomicNum() == 6 and atom.GetIdx() not in visited:
                visited.add(atom.GetIdx())
                chain_atoms.append(atom)
                # Enqueue neighboring carbons if not already visited.
                for nb in atom.GetNeighbors():
                    if nb.GetAtomicNum() == 6 and nb.GetIdx() not in visited:
                        queue.append(nb)
        # Check if this connected carbon subgraph covers all carbons in the molecule.
        if len(chain_atoms) != nC:
            continue
        # For each carbon in the chain (except the aldehyde carbon), check for at least one hydroxyl.
        open_chain_valid = True
        for atom in chain_atoms:
            if atom.GetIdx() == aldehyde_atom.GetIdx():
                continue
            oh_found = False
            for nb in atom.GetNeighbors():
                # Only consider substituents not part of the main chain
                if nb.GetAtomicNum() == 8 and nb.GetIdx() not in [a.GetIdx() for a in chain_atoms]:
                    if is_hydroxyl(nb):
                        oh_found = True
                        break
                # In some cases the –OH might be directly attached to a chain carbon even if that oxygen
                # is also connected to another carbon outside the chain. So also allow oxygens that are not
                # part of the backbone.
                elif nb.GetAtomicNum() == 8:
                    if is_hydroxyl(nb):
                        oh_found = True
                        break
            if not oh_found:
                open_chain_valid = False
                break
        if open_chain_valid:
            return True, "Open-chain aldose pattern detected with terminal aldehyde and hydroxylated carbon chain"
    
    # Attempt 2: Cyclic (hemiacetal) aldose detection
    # Look for 5- or 6-membered rings. Typical cyclic aldoses (furanoses and pyranoses)
    # have one ring oxygen and all ring carbons decorated with an exocyclic hydroxyl.
    ring_info = mol.GetRingInfo().AtomRings()
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        # Identify atoms in the ring by their atomic number.
        oxygens_in_ring = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        carbons_in_ring = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        # For a proper sugar ring, we expect exactly one oxygen.
        if len(oxygens_in_ring) != 1:
            continue
        # For each carbon in the ring, check that it has an exocyclic hydroxyl.
        has_hydroxyl_all = True
        for idx in carbons_in_ring:
            atom = mol.GetAtomWithIdx(idx)
            found_oh = False
            for nb in atom.GetNeighbors():
                # We want substituents that fall outside the ring.
                if nb.GetIdx() in ring:
                    continue
                if nb.GetAtomicNum() == 8 and is_hydroxyl(nb):
                    found_oh = True
                    break
            if not found_oh:
                has_hydroxyl_all = False
                break
        # The number of ring carbons in a sugar ring should be (ring size - 1)
        if has_hydroxyl_all and len(carbons_in_ring) == (len(ring) - 1):
            return True, "Cyclic hemiacetal sugar pattern detected in a ring with proper hydroxyl substitutions"
    
    return False, "No open-chain aldehyde pattern or cyclic hemiacetal sugar pattern detected"

# Example test cases:
if __name__ == "__main__":
    test_cases = [
        # True positives (expected aldoses)
        ("[H][C@](O)(CO)[C@]([H])(O)C=O", "L-erythrose (open-chain)"),
        ("OC[C@@H](O)[C@H](O)[C@H](O)C=O", "aldehydo-D-lyxose (open-chain)"),
        ("O[C@@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O", "beta-D-lyxopyranose (cyclic)"),
        ("OC[C@H]1OC(O)[C@H](O)[C@@H]1O", "D-ribofuranose (cyclic)"),
        ("O1[C@@H]([C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)CO", "beta-D-idopyranose (cyclic)"),
        # False positive examples should return False
        ("O[C@@H]([C@H](O)[C@H](O)C(O)=O)[C@@H](OC)C=O", "2-O-methyl-gluconic acid (should not be aldose)")
    ]
    for smi, name in test_cases:
        result, reason = is_aldose(smi)
        print(f"SMILES: {smi} | {name} -> {result} ; Reason: {reason}")