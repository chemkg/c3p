"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: Indole alkaloid
Definition: An alkaloid containing an indole skeleton.
This improved algorithm searches for a five-membered aromatic ring
(with exactly one nitrogen) that is fused (shares two atoms) with a six-membered aromatic ring.
It also checks that the molecule contains at least one nitrogen atom.
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid is defined as an alkaloid containing an indole skeleton.
    
    This implementation first parses the molecule and then examines its ring system.
    It looks for a five-membered aromatic ring (the pyrrole ring) containing exactly one nitrogen,
    which is fused (sharing exactly two atoms) with a six-membered aromatic ring (benzene ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an indole alkaloid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve all rings (as tuples of atom indices) using RDKit's ring info.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # list of tuples
    
    found_indole = False

    # Iterate over all rings to search for a fused indole system.
    for ring in atom_rings:
        # Look for a possible five-membered ring candidate:
        if len(ring) == 5:
            # Ensure every atom in the ring is aromatic.
            if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                continue
            # Count nitrogen atoms in this five-membered ring.
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_count != 1:
                continue

            # Look for a six-membered aromatic ring that shares exactly two atoms with this candidate ring.
            for ring2 in atom_rings:
                if len(ring2) == 6:
                    if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring2):
                        continue
                    shared_atoms = set(ring) & set(ring2)
                    if len(shared_atoms) == 2:
                        found_indole = True
                        break
            if found_indole:
                break

    if not found_indole:
        return False, "No indole skeleton found"

    # As an additional check, ensure the molecule contains at least one nitrogen atom.
    n_total = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_total == 0:
        return False, "Molecule lacks nitrogen atoms; not an alkaloid"

    return True, "Molecule contains an indole skeleton (a 5-membered aromatic ring fused with a 6-membered aromatic ring) and nitrogen, consistent with an indole alkaloid"

# For testing purposes (if desired), one may include:
if __name__ == "__main__":
    # Example SMILES string for Ochropposinine (a known indole alkaloid)
    test_smiles = "CC[C@]1(CN2CCC=3C4=CC(=C(C=C4NC3[C@@]2(C[C@@]1(CCO)[H])[H])OC)OC)[H]"
    result, reason = is_indole_alkaloid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)