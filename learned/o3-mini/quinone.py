"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: Quinones, defined as compounds having a fully conjugated cyclic dione structure.
For example, benzoquinones and polycyclic/heterocyclic analogues.
This improved heuristic inspects each ring for:
  - A minimum ring size (≥5)
  - A high fraction (≥90%) of ring atoms that are aromatic or sp2‐hybridized
  - Exactly two carbonyl (C(=O)) groups in that ring.
  - Each carbonyl atom has at least 2 in‐ring neighbors that are aromatic/sp2.
  - The two carbonyl atoms are not too “close”, i.e. the shortest bond path, computed
    on the full molecule but restricted to the ring, is at least 3 bonds long.
Note: Given the wide range of quinone structures, this remains only a heuristic.
"""

from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    The heuristic examines each ring in the molecule and searches for a ring (or a fused system)
    having exactly two carbonyl groups (C(=O)) embedded in a conjugated system.
    In addition, for each carbonyl in the ring, at least 2 of the ring neighbors must be sp2‐hybridized or aromatic,
    and the two carbonyls must be separated by at least 3 bonds (within the ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a quinone, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass  # continue even if sanitization fails
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Helper: determine if an atom is a carbonyl carbon (C that has a double bond to O)
    def is_carbonyl(atom):
        if atom.GetAtomicNum() != 6:
            return False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    return True
        return False
    
    # Loop over rings, looking for one that meets the quinone criteria.
    for ring in atom_rings:
        # Ignore small rings
        if len(ring) < 5:
            continue
        
        # Get atoms in the ring.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Check that most of the atoms in the ring are conjugated (sp2 or aromatic)
        n_conjugated = sum(1 for atom in ring_atoms if atom.GetIsAromatic() or atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2)
        if n_conjugated < 0.9 * len(ring_atoms):
            continue
        
        # Identify carbonyl atoms within this ring.
        carbonyl_atoms = [atom for atom in ring_atoms if is_carbonyl(atom)]
        if len(carbonyl_atoms) != 2:
            continue
        
        # For each carbonyl atom, require that at least 2 of its neighbors (that are also in the ring)
        # are sp2 or aromatic to ensure embedded conjugation.
        embedded = True
        for atom in carbonyl_atoms:
            # Get neighbors that are in the ring.
            in_ring_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring]
            n_conj_nbr = sum(1 for nbr in in_ring_neighbors if nbr.GetIsAromatic() or nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP2)
            if n_conj_nbr < 2:
                embedded = False
                break
        if not embedded:
            continue
        
        # Check the separation of the two carbonyls in the ring.
        # Use RDKit's GetShortestPath to compute the bond path between the two carbonyl carbons.
        c1 = carbonyl_atoms[0].GetIdx()
        c2 = carbonyl_atoms[1].GetIdx()
        path = Chem.rdmolops.GetShortestPath(mol, c1, c2)
        # Ensure that the entire path is within the ring.
        if not set(path).issubset(set(ring)):
            # The carbonyls may lie in different rings; skip this ring.
            continue
        # Number of bonds is (number of atoms - 1)
        bond_distance = len(path) - 1
        if bond_distance < 3:
            # Too close (e.g. adjacent carbonyls), likely not a typical quinone.
            continue
            
        # If a ring passes all these criteria, we classify the molecule as a quinone.
        return True, "Found a ring with two carbonyl groups embedded in a fully conjugated cyclic system"
    
    return False, "No fully conjugated cyclic dione (quinone) structure was found"


# For testing purposes (optional)
if __name__ == "__main__":
    # Test a few SMILES strings (both quinones and non-quinones)
    test_smiles = [
        "COC1=C(OC)C(=O)C(C)=CC1=O",  # ubiquinone-0 (expected True)
        "O=C1C=C(OC)[C@](O)(C)C([C@@]1(OC(=O)C(CC)C)C)=O",  # Phomaligol A (expected True)
        "C1(=CC=C2NC(C(C2=C1)=O)=O)Br",  # 5-Bromoisatin (expected False)
    ]
    for s in test_smiles:
        result, reason = is_quinone(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*40}")