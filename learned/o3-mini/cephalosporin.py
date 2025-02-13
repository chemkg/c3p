"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: Cephalosporin
A cephalosporin is characterized by a cephem nucleus: a beta-lactam (4-membered cyclic amide) fused to a 6-membered dihydrothiazine ring (which contains a sulfur).
This program attempts to detect that fused bicyclic system.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    The algorithm looks for a 4-membered beta-lactam ring (one nitrogen, one carbonyl carbon with a double-bonded oxygen)
    fused to a 6-membered ring containing at least one sulfur atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a cephalosporin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use RDKit to retrieve ring information from the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No ring information could be detected"
    
    # Identify candidate beta-lactam rings: 4-membered rings with at least one nitrogen and a carbonyl group.
    beta_lactam_rings = []
    for ring in atom_rings:
        if len(ring) == 4:
            # Retrieve the atoms in this ring.
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Check for at least one nitrogen atom.
            n_nitrogen = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
            # Check for presence of a carbonyl group: a carbon atom in the ring bonded via a double bond to an oxygen.
            has_carbonyl = False
            for atom in atoms_in_ring:
                if atom.GetAtomicNum() == 6:
                    for bond in atom.GetBonds():
                        # Check for double bond to an oxygen.
                        if bond.GetBondTypeAsDouble() == 2.0:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8:
                                has_carbonyl = True
                                break
                    if has_carbonyl:
                        break
            # Accept ring if it contains at least one nitrogen and a carbonyl.
            if n_nitrogen >= 1 and has_carbonyl:
                beta_lactam_rings.append(set(ring))
    
    if not beta_lactam_rings:
        return False, "No beta-lactam (4-membered amide) ring detected"
    
    # Identify candidate 6-membered rings that contain at least one sulfur atom (indicative of the dihydrothiazine ring).
    six_membered_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            if any(atom.GetAtomicNum() == 16 for atom in atoms_in_ring):
                six_membered_rings.append(set(ring))
    
    if not six_membered_rings:
        return False, "No 6-membered ring containing sulfur detected (dihydrothiazine ring missing)"
    
    # Check for fused rings: the beta-lactam ring and the 6-membered ring must share at least 2 atoms.
    for beta_ring in beta_lactam_rings:
        for six_ring in six_membered_rings:
            common_atoms = beta_ring.intersection(six_ring)
            if len(common_atoms) >= 2:
                return True, "Molecule has a beta-lactam ring fused with a 6-membered sulfur-containing ring (cephem nucleus)"
    
    return False, "No fused beta-lactam and 6-membered sulfur-containing ring system detected"

# Example usage:
if __name__ == '__main__':
    # An example SMILES from the provided list of cephalosporins:
    test_smiles = "CC1=C(C(=NN1CC(C)C(=O)NC2C3N(C2=O)C(=C(CS3)COC(=O)C)C(=O)O)C)Cl"
    result, reason = is_cephalosporin(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)