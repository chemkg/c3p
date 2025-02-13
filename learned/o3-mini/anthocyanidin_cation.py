"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: Anthocyanidin Cation (flavylium‐type aglycon of anthocyanin cation)
Definition: Any organic cation that is an aglycon of anthocyanin cation; they are oxygenated derivatives of flavylium (2-phenylchromenylium).
This heuristic function checks for a positively charged oxygen within a six-membered (pyran) ring that is fused to at least one other aromatic ring,
and an aromatic (phenyl) substituent attached to one of the ring carbons.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation (flavylium-type aglycon) based on its SMILES string.
    
    Heuristic criteria:
      1. The molecule must be valid.
      2. The overall molecule must carry positive charge(s).
      3. There must be at least one aromatic oxygen with formal charge +1.
      4. This oxygen should lie within a six-membered aromatic ring (the pyran ring of the flavylium core).
      5. The ring containing the [O+ ] should be fused to at least one other aromatic ring.
      6. At least one carbon of the candidate ring should have a substituent that is an aromatic ring (a phenyl group)
         representing the 2-phenyl substitution expected in a 2-phenylchromenylium core.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an anthocyanidin cation, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure explicit hydrogens are added to support aromaticity assignments if needed.
    mol = Chem.AddHs(mol)
    Chem.SanitizeMol(mol)
    
    # Check that the molecule has at least one positive formal charge.
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge <= 0:
        return False, "Molecule does not carry a net positive charge"
    
    # Look for candidate atoms: aromatic oxygen atoms with formal charge +1.
    candidate_oatoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetFormalCharge() == 1 and atom.GetIsAromatic():
            candidate_oatoms.append(atom.GetIdx())
    
    if not candidate_oatoms:
        return False, "No aromatic oxygen atom with +1 charge found (flavylium moiety missing)"
    
    # Get ring information.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Function to check if a given ring (given as a tuple of atom indices) is aromatic.
    def is_ring_aromatic(ring):
        return all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
    
    # Check candidates: look for a six-membered aromatic ring that contains a candidate [O+] atom.
    candidate_rings = []
    for ring in ring_info:
        if len(ring) == 6 and is_ring_aromatic(ring):
            if any(idx in ring for idx in candidate_oatoms):
                candidate_rings.append(ring)
    
    if not candidate_rings:
        return False, "No six-membered aromatic ring containing an [O+ ] atom found"
    
    # For each candidate ring, check for ring fusion: at least one other aromatic ring sharing 2 atoms.
    fused_found = False
    for ring in candidate_rings:
        for other_ring in ring_info:
            if ring == other_ring:
                continue
            # Check if rings share at least 2 atoms (typical for fused rings)
            shared_atoms = set(ring).intersection(set(other_ring))
            if len(shared_atoms) >= 2 and is_ring_aromatic(other_ring):
                fused_found = True
                break
        if fused_found:
            candidate_ring = ring
            break
    if not fused_found:
        return False, "No fused aromatic ring detected with the candidate flavylium ring"
    
    # Check for a phenyl substituent (an aromatic ring attached as a side chain)
    # We look for a carbon atom in the candidate ring that has a neighbor (not in candidate_ring) that is aromatic,
    # and that neighbor is part of a six-membered ring (phenyl group).
    phenyl_found = False
    candidate_ring_set = set(candidate_ring)
    for idx in candidate_ring:
        atom = mol.GetAtomWithIdx(idx)
        # Consider non-oxygen atoms in the candidate ring (should be sp2 carbons).
        if atom.GetSymbol() != 'C':
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in candidate_ring_set:
                continue
            if nbr.GetIsAromatic() and nbr.GetSymbol() == 'C':
                # Check if the neighbor is part of a six-membered aromatic ring.
                for ring in ring_info:
                    if len(ring) == 6 and nbr.GetIdx() in ring and is_ring_aromatic(ring):
                        phenyl_found = True
                        break
            if phenyl_found:
                break
        if phenyl_found:
            break
    
    if not phenyl_found:
        return False, "No phenyl substituent attached to the flavylium core detected"
    
    # If all heuristic criteria are met, we consider the molecule as an anthocyanidin cation.
    return True, "Molecule contains a flavylium-type (2-phenylchromenylium) core with a positively charged oxygen, fused aromatic rings, and an attached phenyl group"

# For complex cases that may not be unambiguously classified, one could instead return (None, None).
# For example, if you prefer to use a more conservative approach, uncomment the next line.
# return None, None

# Example usage:
if __name__ == '__main__':
    # You can test one of the provided examples (e.g., hirsutidin)
    test_smiles = "COc1cc(O)c2cc(O)c([o+]c2c1)-c1cc(OC)c(O)c(OC)c1"
    result, reason = is_anthocyanidin_cation(test_smiles)
    print("Result:", result)
    print("Reason:", reason)