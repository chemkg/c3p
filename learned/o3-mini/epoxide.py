"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: Epoxide – Any cyclic ether in which the oxygen atom forms part of a 3‐membered ring.
For our purposes we only classify an epoxide if the 3‐membered cyclic ether is an isolated ring,
i.e. none of the ring atoms participate in any other ring.
"""

from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is defined as any cyclic ether in which the oxygen atom forms part of a 
    3-membered ring that is isolated (i.e. not fused to other rings).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an epoxide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # list of tuples (each tuple is indices of atoms in a ring)

    if not atom_rings:
        return False, "No rings found in molecule"
    
    # Precompute for each atom how many rings it is in.
    atom_ring_count = {atom.GetIdx(): 0 for atom in mol.GetAtoms()}
    for ring in atom_rings:
        for idx in ring:
            atom_ring_count[idx] += 1

    # Loop through each ring and check if any ring is a 3-membered cyclic ether
    for ring in atom_rings:
        if len(ring) == 3:
            # Get atomic numbers for atoms in the ring
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            atomic_nums = sorted([atom.GetAtomicNum() for atom in atoms])
            # For an epoxide, we expect two carbons (atomic number 6) and one oxygen (atomic number 8)
            if atomic_nums == [6, 6, 8]:
                # Check if this ring is isolated (none of its atoms lie in another ring)
                if all(atom_ring_count[idx] == 1 for idx in ring):
                    return True, "Found an isolated 3-membered cyclic ether (epoxide) group"
    
    # If no suitable isolated 3-membered ring is found, then it fails the epoxide criteria.
    return False, "Does not contain an isolated 3-membered cyclic ether (epoxide) group"

# Example usage (you can remove or comment out these lines when integrating into a larger codebase)
if __name__ == "__main__":
    examples = [
        ("ClCC1CO1", "epichlorohydrin (should be True)"),
        ("C1CC1", "Cyclopropane (should be False)"),
        ("O1C(C)C1", "Epoxide-like ring with extra substituent (likely True if isolated)")
    ]
    for smi, desc in examples:
        result, reason = is_epoxide(smi)
        print(f"SMILES: {smi} | {desc} --> {result}, {reason}")