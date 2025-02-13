"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: Cephalosporin
A cephalosporin is defined as having a cephem nucleus – a beta‐lactam (4‐membered
cyclic amide) fused to a 6‐membered heterocycle (typically a dihydrothiazine ring, but in some cases
an oxygen may replace sulfur). This implementation attempts to detect that bicyclic core after cleaning
the molecule from salts or small fragments.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    
    This function performs the following steps:
      1. Parses the SMILES and removes salts/water by retaining the largest fragment.
      2. Analyzes the rings in the molecule.
      3. Identifies candidate beta-lactam rings – 4-member rings that have at least one nitrogen and at least one carbonyl group.
      4. Identifies candidate 6-membered rings – rings of length 6 that contain at least one heteroatom (S or O).
      5. Checks for a fused system by ensuring that one 4-membered ring and one 6-membered ring share at least 2 atoms.
      
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
    
    # Retrieve ring information for the cleaned molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No ring systems detected in the molecule"
    
    # Identify candidate beta-lactam rings: 4-membered rings with at least one nitrogen and a carbonyl group.
    beta_lactam_rings = []
    for ring in atom_rings:
        if len(ring) == 4:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            n_nitrogen = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
            has_carbonyl = False
            # Look for a carbon (atomic number 6) in the ring that is double-bonded to an oxygen.
            for atom in atoms_in_ring:
                if atom.GetAtomicNum() == 6:
                    for bond in atom.GetBonds():
                        # bond.GetBondTypeAsDouble() returns 2.0 for double bonds.
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
    
    # Identify candidate 6-membered rings containing at least one heteroatom.
    # Typically cephalosporins have a dihydrothiazine (contains S) but some variants have O.
    six_membered_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Check if any atom in ring is heteroatom (S (16) or O (8)) and not carbon.
            if any(atom.GetAtomicNum() in (8, 16) for atom in atoms_in_ring):
                six_membered_rings.append(set(ring))
    
    if not six_membered_rings:
        return False, "No 6-membered heterocycle (dihydrothiazine or related ring) detected"
    
    # Check for fused rings: the beta-lactam ring must be fused (share at least 2 atoms) with a 6-membered ring.
    for beta_ring in beta_lactam_rings:
        for six_ring in six_membered_rings:
            common_atoms = beta_ring.intersection(six_ring)
            if len(common_atoms) >= 2:
                return True, "Molecule has a beta-lactam ring fused with a 6-membered heterocycle (cephem nucleus)"
    
    return False, "No fused system of a beta-lactam ring and a 6-membered heterocycle was detected"


# Example usage:
if __name__ == '__main__':
    # Test with one of the provided cephalosporin SMILES strings, e.g. cefpodoxime.
    test_smiles = "[H][C@]12SCC(COC)=C(N1C(=O)[C@H]2NC(=O)C(=N/OC)\\c1csc(N)n1)C(O)=O"
    result, reason = is_cephalosporin(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)