"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: Organometalloidal Compound
Definition: A compound having bonds between one or more metalloid (here we focus on arsenic, As)
and one or more carbon atoms belonging to an organyl group.
Heuristic:
  – First, we find arsenic atoms (atomic number 33).
  – For each As atom, we get its carbon neighbors.
  – If only one carbon is attached, we accept (e.g. methylarsonic acid).
  – If more than one carbon is attached, we require that at least one carbon is not “simple methyl.”
    Here, a carbon is considered “simple methyl” if it is aliphatic (non‐aromatic) and (aside from the As)
    has only one heavy-atom neighbor.
  – For carbons that are aromatic, we further check that the carbon belongs to a “simple” benzene ring,
    i.e. (a) the ring is exactly six members, (b) all atoms in the ring are carbons, and (c) the carbon
    does not belong to more than one ring (which may indicate a fused system).
If none of the As–C bonds meets these criteria, the compound is not classified as an organometalloidal compound.
Note: This is only a heuristic.
"""

from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    We focus on arsenic-containing compounds. The molecule is classified as an organometalloidal compound
    if it contains at least one As–C bond where, in the case of multiple carbon substituents, at least one
    carbon is not a simple methyl group. For aromatic carbons the substituent must belong to a nonfused,
    benzene-like ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule qualifies, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve ring information from the molecule
    ring_info = mol.GetRingInfo().AtomRings()
    
    AS_ATOMIC_NUM = 33  # atomic number for arsenic
    
    def is_simple_methyl(carbon, attached_to):
        """
        Checks if the given carbon atom is a "simple methyl" group.
        A carbon is considered a methyl if:
          - It is aliphatic (non-aromatic).
          - Excluding the connection to the As atom, it has exactly one heavy-atom neighbor.
        """
        if carbon.GetSymbol() != "C":
            return False
        if carbon.GetIsAromatic():
            return False  # aromatic carbons are considered more complex
        # Count heavy (non-hydrogen) neighbors excluding the attached metal atom.
        neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != attached_to.GetIdx()]
        return len(neighbors) == 1

    def qualifies_aromatic(carbon):
        """
        For aromatic carbons, we require that the atom is part of a single six-membered ring
        in which every member is a carbon.
        This is our rough way to capture a simple phenyl (organyl) group.
        """
        # Find all rings that contain this atom.
        rings_here = [r for r in ring_info if carbon.GetIdx() in r]
        # Reject if the atom appears in more than one ring (i.e. likely in a fused system).
        if len(rings_here) != 1:
            return False
        ring = rings_here[0]
        if len(ring) != 6:
            return False
        # Check that every atom in the ring is a carbon.
        for idx in ring:
            if mol.GetAtomWithIdx(idx).GetSymbol() != "C":
                return False
        return True

    # Loop over all atoms looking for arsenic
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != AS_ATOMIC_NUM:
            continue
        
        # Get carbon neighbors for this arsenic atom.
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not carbon_neighbors:
            continue  # this As atom is not bonded to any carbon
        
        # If there is exactly one carbon neighbor, we accept the bond (e.g. methylarsonic acid).
        if len(carbon_neighbors) == 1:
            neigh = carbon_neighbors[0]
            return True, f"Found a bond between As and C ({neigh.GetSymbol()}). Only one carbon substituent present."
        
        # More than one carbon is bonded.
        # We now require that at least one of these carbons is not a simple methyl.
        for cn in carbon_neighbors:
            if cn.GetIsAromatic():
                # For aromatic carbon, require that it belongs to a simple benzene ring.
                if qualifies_aromatic(cn):
                    return True, f"Found a bond between As and aromatic C from a benzene-like ring."
            else:
                # For aliphatic carbons, check if it is more than a methyl.
                if not is_simple_methyl(cn, attached_to=atom):
                    return True, f"Found a bond between As and an aliphatic C that is not a simple methyl."
        # If all attached carbons are simple methyl groups, then we do not qualify.
        return False, "Arsenic is bonded to multiple only-methyl groups; does not qualify as an organometalloidal compound."
    
    return False, "No valid As–C bond (with proper organyl group) detected."

# Example usage (testing some of the examples):
if __name__ == "__main__":
    test_examples = [
        # True positives
        ("Nc1cc(ccc1O)[As]1[As]([As]1c1ccc(O)c(N)c1)c1ccc(O)c(N)c1", "arsphenamine trimer"),
        ("[As](=O)(CCCCCCCCCCCCCCCCC)(C)C", "1-dimethylarsinoyl-heptadecane"),
        ("OC(=O)C[As](O)(O)=O", "arsenoacetic acid"),
        ("C[As](O)(O)=O", "methylarsonic acid"),
        ("O[As](O)(=O)c1ccccc1", "phenylarsonic acid"),
        # Cases that should not qualify (e.g. only methyl groups)
        ("C[As](C)(O)=O", "dimethylarsinic acid"),
        ("C[As](C)O", "dimethylarsinous acid"),
        # One false positive example from previous attempt (ideally now rejected):
        ("Oc1ccc2c(oc3c([As]4SCCS4)c(O)ccc3c2=O)c1[As]1SCCS1", "HOxAsH-bis(1,2-ethanedithiol)")
    ]
    
    for smi, name in test_examples:
        res, reason = is_organometalloidal_compound(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nRESULT: {res}\nREASON: {reason}\n{'-'*60}")