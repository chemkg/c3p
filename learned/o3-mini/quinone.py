"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: Quinones, defined as compounds having a fully conjugated cyclic dione structure,
for example benzoquinones (polycyclic and heterocyclic analogues are included).
"""

from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    The function looks for a cyclic structure (ring) that contains exactly two carbonyl groups (C(=O))
    and is largely conjugated (most ring atoms are sp2 or aromatic).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a quinone, False otherwise.
        str: The reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure Kekulization (so that double bonds are recognized correctly)
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass  # If kekulization fails, we may still be able to work with the molecule

    # Retrieve all the rings in the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Helper: is a given atom a carbonyl carbon? It must be carbon and have a double bonded oxygen.
    def is_carbonyl(atom):
        if atom.GetAtomicNum() != 6:
            return False
        for neighbor in atom.GetNeighbors():
            # Check that neighbor is oxygen
            if neighbor.GetAtomicNum() == 8:
                # Look for double bond connecting atom and neighbor
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    return True
        return False
    
    # Iterate over each ring in the molecule
    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count how many carbonyl carbons are in this ring
        carbonyl_count = sum(1 for atom in ring_atoms if is_carbonyl(atom))
        
        # We are looking for a cyclic dione structure
        if carbonyl_count == 2:
            # Next, check if the ring is "fully conjugated".
            # A simple heuristic: count atoms that are either aromatic or have sp2 hybridization.
            sp2_or_aromatic = 0
            for atom in ring_atoms:
                # Using both aromatic flag and hybridization as clues.
                if atom.GetIsAromatic() or atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                    sp2_or_aromatic += 1
            # If most ring atoms are conjugated, consider the ring as fully conjugated.
            if sp2_or_aromatic >= 0.8 * len(ring_atoms):
                return True, "Found a ring with two carbonyl groups in a conjugated cyclic system"
    
    return False, "No fully conjugated cyclic dione (quinone) structure was found"