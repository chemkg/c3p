"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: Cephalosporin

A cephalosporin is defined as having a cephem nucleus – a beta‐lactam (4‐membered cyclic amide)
fused to a 6‐membered heterocycle (typically a dihydrothiazine ring, sometimes with oxygen instead of sulfur).
This implementation attempts to detect that bicyclic core after cleaning the molecule from salts or small fragments.
An important improvement compared to our previous version is that instead of accepting any fusion (≥2 common atoms)
between a 4-membered ring and a 6-membered ring, we require that the two rings share exactly 2 atoms,
producing an 8-atom fused system. This is more in line with the cephem nucleus.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    
    Steps:
      1. Parse the SMILES and remove salts/solvents by retaining the largest fragment.
      2. Obtain ring information.
      3. Identify candidate beta-lactam rings: 4-membered rings that contain at least one nitrogen and a carbonyl group.
      4. Identify candidate six-membered rings: 6-membered rings that contain at least one heteroatom (S or O).
      5. For each candidate pair, require that their union has exactly 8 atoms (i.e. they share exactly 2 atoms),
         which corresponds to the cephem nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a cephalosporin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove salts/solvents by retaining only the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments could be extracted from molecule"
    mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())

    # Retrieve ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No ring systems detected in the molecule"

    # Identify candidate beta-lactam rings: 4-membered rings with >=1 nitrogen and a carbonyl group.
    beta_lactam_rings = []
    for ring in atom_rings:
        if len(ring) == 4:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            n_nitrogen = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
            has_carbonyl = False
            # Check for a carbon atom in the ring having a double bond to an oxygen (carbonyl)
            for atom in atoms_in_ring:
                if atom.GetAtomicNum() == 6:
                    for bond in atom.GetBonds():
                        # bond.GetBondTypeAsDouble() returns 2.0 for a double bond.
                        if bond.GetBondTypeAsDouble() == 2.0:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8:
                                has_carbonyl = True
                                break
                    if has_carbonyl:
                        break
            if n_nitrogen >= 1 and has_carbonyl:
                beta_lactam_rings.append(set(ring))

    if not beta_lactam_rings:
        return False, "No beta-lactam (4-membered cyclic amide) ring detected"

    # Identify candidate 6-membered rings that contain at least one heteroatom (S or O).
    six_membered_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Check if any atom in ring is sulfur (16) or oxygen (8)
            if any(atom.GetAtomicNum() in (8, 16) for atom in atoms_in_ring):
                six_membered_rings.append(set(ring))

    if not six_membered_rings:
        return False, "No 6-membered heterocycle (dihydrothiazine or related ring) detected"

    # Check for correctly-fused rings.
    # In a cephem nucleus, the 4-member ring and the 6-member ring share exactly 2 atoms.
    # (i.e. the union is 4 + 6 – 2 = 8 atoms)
    for beta_ring in beta_lactam_rings:
        for six_ring in six_membered_rings:
            common_atoms = beta_ring.intersection(six_ring)
            union_atoms = beta_ring.union(six_ring)
            if len(common_atoms) == 2 and len(union_atoms) == 8:
                return True, "Molecule has a beta-lactam ring fused with a 6-membered heterocycle (cephem nucleus)"
    
    return False, "No fused system with an 8-atom cephem nucleus was detected"


# Example usage:
if __name__ == '__main__':
    # Test with one of the provided cephalosporin SMILES strings, e.g. cefpodoxime.
    test_smiles = "[H][C@]12SCC(COC)=C(N1C(=O)[C@H]2NC(=O)C(=N/OC)\\c1csc(N)n1)C(O)=O"
    result, reason = is_cephalosporin(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)