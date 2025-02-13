"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: Guaiacols
Defined as: Any phenol carrying an additional methoxy substituent at the ortho‐position.
This program checks each aromatic benzene ring in the molecule and looks for an adjacent pair of substituents:
one being a hydroxyl group (-OH) and the other a methoxy group (-OCH3).
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
        str: Explanation of the classification result.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens so we can reliably count hydrogens (both implicit and explicit)
    mol = Chem.AddHs(mol)
    
    # Helper function to check if an atom (outside the ring) is a hydroxyl group
    def has_hydroxy(atom, ring_set):
        # Look at substituents of the ring atom that are NOT in the ring
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in ring_set:
                continue
            # Check if the neighbor is an oxygen.
            if nbr.GetSymbol() == "O":
                # For a hydroxyl, oxygen must have at least one hydrogen.
                # Use GetTotalNumHs() to count implicit+explicit hydrogens.
                if nbr.GetTotalNumHs() >= 1:
                    return True
        return False
    
    # Helper function to check if an atom (outside the ring) is a methoxy group (-OCH3)
    def has_methoxy(atom, ring_set):
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in ring_set:
                continue
            if nbr.GetSymbol() == "O":
                # A methoxy oxygen should be connected to exactly one non‐ring carbon.
                non_ring_nbs = [n for n in nbr.GetNeighbors() if n.GetIdx() not in ring_set]
                if len(non_ring_nbs) != 1:
                    continue
                carbon = non_ring_nbs[0]
                if carbon.GetAtomicNum() != 6:
                    continue
                # Check that the carbon (in the methoxy group) has exactly 3 hydrogens (total, not only explicit)
                # This typically means it is a methyl group.
                if carbon.GetTotalNumHs() == 3:
                    return True
        return False

    # Get ring information from the molecule and loop over the rings.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # We only consider six-membered rings (benzene-like) that are aromatic.
        if len(ring) != 6:
            continue
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        ring_set = set(ring)
        # Iterate over each pair of adjacent atoms in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Loop over adjacent ring atoms.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in ring_set:
                    continue
                # To avoid checking the same pair twice, enforce an order.
                if idx < nbr.GetIdx():
                    # Check if one atom has -OH and the other has -OCH3 (ortho substituents).
                    if (has_hydroxy(atom, ring_set) and has_methoxy(nbr, ring_set)) or \
                       (has_methoxy(atom, ring_set) and has_hydroxy(nbr, ring_set)):
                        return True, "Contains a guaiacol moiety (phenol with an ortho-methoxy substituent)"
    
    return False, "Does not contain a guaiacol moiety: missing phenol with ortho-methoxy substituent"

# Example usage (if run as main, only for testing purposes)
if __name__ == "__main__":
    # A few test SMILES: guaiacol itself and one that is not.
    test_molecules = [
        ("COc1ccccc1O", "guaiacol"),
        ("c1ccccc1O", "phenol without methoxy")
    ]
    for smi, name in test_molecules:
        result, reason = is_guaiacols(smi)
        print(f"SMILES: {smi}  NAME: {name}  -> {result}  REASON: {reason}")