"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: Organometalloidal Compound
Definition: A compound having bonds between one or more metalloid atoms (here we focus on arsenic, As)
and one or more carbon atoms belonging to an organyl group.
Heuristic:
  – Find all arsenic atoms.
  – For each As, get its carbon neighbors.
  – If exactly one carbon is attached, we accept (e.g. methylarsonic acid).
  – If more than one carbon is attached, at least one must be more than a simple methyl.
       For aliphatic carbons we “measure” the chain length (a true organyl group should be long enough, here ≥16).
       For aromatic carbons we require that the carbon is in a “simple” ring (either 5‐ or 6–membered, and non‐fused).
If none of the As–C bonds meets these criteria, the compound is not classified as organometalloidal.
Note: This heuristic is imperfect.
"""
from rdkit import Chem

# Helper function to decide if an aliphatic carbon is a simple methyl group.
def is_simple_aliphatic(carbon, attached_to):
    """
    Returns True if the (non-aromatic) carbon is only a CH3 – i.e. apart from its bond to the metal,
    it is not bonded to any other heavy (atomic number > 1) atoms.
    """
    if carbon.GetSymbol() != "C":
        return False
    if carbon.GetIsAromatic():
        return False  # don’t check aromatic here
    # Count heavy-atom neighbors other than the attached metalloid.
    count = 0
    for nbr in carbon.GetNeighbors():
        if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != attached_to.GetIdx():
            count += 1
    # A methyl carbon (CH3) has no heavy atoms besides the one (As) it’s attached to.
    return count == 0

# Helper to compute the longest chain length (number of carbon atoms) 
# starting from a given aliphatic carbon, following only non-aromatic carbons.
def get_max_aliphatic_chain_length(atom, coming_from, visited):
    # Only traverse if the carbon is non-aromatic and is truly aliphatic.
    visited.add(atom.GetIdx())
    max_length = 1  # count this atom
    for nbr in atom.GetNeighbors():
        # Avoid going back to the atom we came from or revisiting atoms.
        if nbr.GetIdx() == coming_from.GetIdx() or nbr.GetIdx() in visited:
            continue
        if nbr.GetSymbol() == "C" and not nbr.GetIsAromatic():
            # Proceed recursively.
            branch_length = 1 + get_max_aliphatic_chain_length(nbr, atom, visited.copy())
            if branch_length > max_length:
                max_length = branch_length
    return max_length

# Helper to check if an aromatic carbon is in a “simple” ring.
def qualifies_aromatic(carbon, mol):
    ring_info = mol.GetRingInfo().AtomRings()
    # Get rings that contain this carbon.
    rings_here = [r for r in ring_info if carbon.GetIdx() in r]
    if not rings_here:
        return False
    # For acceptance we require at least one ring that is either five- or six-membered,
    # and that the carbon appears in only that ring (i.e. non-fused).
    for r in rings_here:
        if len(r) in {5, 6}:
            # Check that all atoms in the ring are carbons.
            if all(mol.GetAtomWithIdx(idx).GetSymbol() == "C" for idx in r):
                # Check that carbon is not in more than one ring.
                count = sum(1 for ring in rings_here if carbon.GetIdx() in ring)
                if count == 1:
                    return True
    return False

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    Focus is on arsenic-containing compounds. The molecule is classified as an organometalloidal compound
    if it contains at least one As–C bond for which:
      – if only one C is present, it is accepted (e.g. methylarsonic acid);
      – if multiple C’s are present, at least one C is “non-simple”, defined as:
            * For aliphatic C: the longest carbon chain (starting at the attached C, following only non-aromatic C)
              has length >= 16.
            * For aromatic C: the C is part of a non-fused simple ring (of size 5 or 6).
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule qualifies, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    AS_ATOMIC_NUM = 33  # arsenic

    # Loop over atoms to find arsenic atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != AS_ATOMIC_NUM:
            continue
        # Get all carbon neighbors attached to this As.
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not carbon_neighbors:
            continue  # no As–C bond here
        
        # If only one C substituent is attached, accept immediately.
        if len(carbon_neighbors) == 1:
            cn = carbon_neighbors[0]
            return True, f"Found a bond between As and C ({cn.GetSymbol()}). Only one carbon substituent present."
        
        # More than one carbon substituent: require at least one that is not “simple”.
        for cn in carbon_neighbors:
            # Check aliphatic carbons first.
            if not cn.GetIsAromatic():
                # If it is a simple methyl group, skip.
                if is_simple_aliphatic(cn, attached_to=atom):
                    continue
                # Else, compute the maximum chain length starting from this carbon.
                chain_length = get_max_aliphatic_chain_length(cn, coming_from=atom, visited=set())
                # For our heuristic, we require a sufficiently long chain to be considered an organyl group.
                if chain_length >= 16:
                    return True, f"Found a bond between As and an aliphatic C chain of length {chain_length}."
            else:
                # For aromatic carbons, require that it is part of a non-fused simple ring.
                if qualifies_aromatic(cn, mol):
                    return True, "Found a bond between As and an aromatic C from a simple ring."
        # If none of the attached carbons passes our tests for a non-simple substituent,
        # then this arsenic atom fails the criteria.
        return False, "Arsenic is bonded only to simple (methyl-like) groups or groups that are too small."
    
    return False, "No valid As–C bond (with proper organyl group) detected."

# Example usage (only for testing; in production you may remove or modify these as needed):
if __name__ == "__main__":
    test_examples = [
        # True positives:
        ("Nc1cc(ccc1O)[As]1[As]([As]1c1ccc(O)c(N)c1)c1ccc(O)c(N)c1", "arsphenamine trimer"),
        ("[As](=O)(CCCCCCCCCCCCCCCCC)(C)C", "1-dimethylarsinoyl-heptadecane"),
        ("[As](=O)(CCCCCCCCCCCCCCCCCCCCCCC)(C)C", "1-dimethylarsinoyl-tricosane"),
        ("CC(=O)N[C@@H](Cc1ccc(O)c(c1)\\N=N\\c1ccc(cc1)[As](O)(O)=O)C(=O)NCC(=O)NCC(=O)NNC(=O)OC(C)(C)C", "3-[(4-arsonophenyl)diazenyl]-AcTyrGlyGlyNHNHBoc"),
        ("C[As](O)(O)=O", "methylarsonic acid"),
        ("O[As](O)(=O)c1ccccc1", "phenylarsonic acid"),
        ("Nc1cc(ccc1O)[As]1[As]([As]([As]([As]1c1ccc(O)c(N)c1)c1ccc(O)c(N)c1)c1ccc(O)c(N)c1)c1ccc(O)c(N)c1", "arsphenamine pentamer"),
        ("OC(=O)C[As](O)(O)=O", "arsenoacetic acid"),
        ("C[As]([O-])([O-])=O", "methylarsonate(2-)"),
        ("O[As](O)(=O)c1ccc(cc1)[N+]([O-])=O", "nitarsone"),
        ("C[As](O)([O-])=O", "methylarsonate(1-)"),
        ("NCC[As](O)(O)=O", "2-Aminoethylarsonate"),
        ("[As](=O)(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)(C)C", "1-dimethylarsinoyl-(3Z,6Z,9Z,12Z,15Z,18Z)-docosahexaene"),
        ("[As+](CCO)(C)(C)C", "arsenocholine"),
        ("CC(=O)N[C@@H](Cc1ccc(O)c(c1)\\N=N\\c1ccc(cc1)[As](O)(O)=O)C(=O)NCC(=O)NCC(O)=O", "3-[(4-arsonophenyl)diazenyl]-AcTyrGlyGly"),
        ("Nc1ccc(cc1)[As](O)(O)=O", "arsanilic acid"),
        ("Oc1ccc(cc1[N+]([O-])=O)[As](O)(O)=O", "roxarsone"),
        ("OC(=O)C\\[As]=[As]\\CC(O)=O", "arsenoacetic acid"),
        ("O[As](O)(=O)c1ccc(cc1)[N+]#N", "ABA diazonium"),
        ("[Na+].Nc1ccc(cc1)[As](O)([O-])=O", "sodium arsanilate"),
        ("C[As](C)(O)=O", "dimethylarsinic acid"),
        ("O[As](O)(=O)c1ccc(cc1)\\N=N\\c1ccc(cc1)[As](O)(O)=O", "4,4'-azodibenzenearsonic acid"),
        ("C[As](C)O", "dimethylarsinous acid"),
        ("C[As](C)c1ccccc1[As](C)C", "1,2-phenylenebis(dimethylarsane)"),
        ("O[As]([O-])(=O)c1ccccc1", "phenylarsonate(1-)"),
        ("C[As](C)(C)(C)C", "pentamethyl-lambda(5)-arsane"),  # likely false but included for testing
        # False positives (should be rejected):
        ("Oc1ccc2c(oc3c([As]4SCCS4)c(O)ccc3c2=O)c1[As]1SCCS1", "HOxAsH-bis(1,2-ethanedithiol)"),
        ("[As](=O)(CCCCCCCCCCCCCCC(O)=O)(C)C", "15-dimethylarsinoyl pentadecanoic acid"),
        ("N[C@@H](CCC(=O)N[C@@H](CSCC(=O)Nc1ccc(cc1)[As](O)O)C(=O)NCC(O)=O)C(O)=O", "S-{2-[4-(dihydroxyarsino)phenylamino]-2-oxoethyl}-glutathione"),
        ("C1C=C[As]=C1", "3H-arsole"),
        ("C1CC[AsH]C1", "arsolane"),
        ("C1C=CC=[As]1", "2H-arsole"),
        ("C[As](C)C", "trimethylarsine")
    ]
    
    for smi, name in test_examples:
        res, reason = is_organometalloidal_compound(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nRESULT: {res}\nREASON: {reason}\n{'-'*60}")