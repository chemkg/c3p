"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: Guaiacols
Defined as: Any phenol carrying an additional methoxy substituent at the ortho‐position.
The program scans each aromatic benzene ring in the molecule and checks for an adjacent pair of substituents,
with one being a hydroxyl group (-OH) and the other a methoxy group (-OCH3).
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is defined as any phenol molecule carrying an additional methoxy group (–OCH3)
    on an ortho-position relative to the –OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a guaiacol moiety, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens so we can check for -OH hydrogens and CH3 groups
    mol = Chem.AddHs(mol)
    
    # Helper function to check if a ring atom carries an -OH group.
    def has_hydroxy(atom, ring_set):
        # Loop over neighbors that are not part of the ring.
        for nb in atom.GetNeighbors():
            if nb.GetIdx() in ring_set:
                continue
            # Check if the neighbor is an oxygen.
            if nb.GetAtomicNum() == 8:
                # Check if this oxygen has at least one hydrogen attached.
                for nnb in nb.GetNeighbors():
                    if nnb.GetAtomicNum() == 1:
                        return True
        return False
    
    # Helper function to check if a ring atom carries an -OCH3 group.
    def has_methoxy(atom, ring_set):
        for nb in atom.GetNeighbors():
            if nb.GetIdx() in ring_set:
                continue
            # Look for an oxygen substituent.
            if nb.GetAtomicNum() == 8:
                # For a methoxy, the oxygen should be attached to the ring atom and exactly one carbon.
                # Check that oxygen has exactly two neighbors: one must be the ring carbon and the other a carbon from CH3.
                if nb.GetDegree() != 2:
                    continue
                # Identify the neighbor that is not in the ring.
                non_ring_nbs = [n for n in nb.GetNeighbors() if n.GetIdx() not in ring_set]
                if len(non_ring_nbs) != 1:
                    continue
                carbon = non_ring_nbs[0]
                if carbon.GetAtomicNum() != 6:
                    continue
                # Check that the carbon is CH3. (After AddHs, explicit H count should be 3.)
                # Some molecules may show fewer explicit hydrogens if using implicit hydrogens, so we use GetNumExplicitHs.
                if carbon.GetNumExplicitHs() == 3:
                    return True
        return False

    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    # Loop over all rings; we only consider rings that are benzene-like (6 atoms, aromatic)
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        # Check that all atoms in the ring are aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        ring_set = set(ring)
        # For each bond (adjacent pair) in the ring, check if one atom has -OH and its neighbor has -OCH3.
        # We consider each pair only once.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check neighbors that are part of the ring.
            for nb in atom.GetNeighbors():
                if nb.GetIdx() not in ring_set:
                    continue
                # To avoid double checking, only consider pairs where our index is lower.
                if idx < nb.GetIdx():
                    # Check for the ortho pattern: one has hydroxyl and the other has methoxy.
                    if (has_hydroxy(atom, ring_set) and has_methoxy(nb, ring_set)) or (has_methoxy(atom, ring_set) and has_hydroxy(nb, ring_set)):
                        return True, "Contains a guaiacol moiety (phenol with an ortho-methoxy substituent)"
    
    return False, "Does not contain a guaiacol moiety: missing phenol with ortho-methoxy substituent"