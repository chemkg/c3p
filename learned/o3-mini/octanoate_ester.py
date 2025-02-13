"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: Octanoate ester
Definition: Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).

The algorithm:
  1. Parse the SMILES.
  2. Look for ester groups using the SMARTS pattern "[#6X3](=O)[OX2H0]".
  3. For each ester found:
       a. Identify the acyl branch (the carbon attached to the carbonyl aside from the ester oxygen).
       b. “Walk” the acyl chain following only single bonds between sp3 carbons.
       c. Check that the chain has exactly eight carbons (including the carbonyl carbon) and that its terminal carbon is a methyl group.
       d. Check that the ester oxygen (on the alcohol side) is not connected (even within one bond) to a phosphorus atom (to help filter out e.g. phospholipids).
       e. Compute the “relative size” of the putative octanoate fragment (counting the atoms in the linear chain plus the carbonyl oxygen) compared with the heavy atoms in the entire molecule. If the octanoate accounts for at least ~35% of the heavy atoms, then we classify the molecule as an octanoate ester.
  4. If any ester match passes these criteria, return True.
  
Note: This is a heuristic approach and may still be imperfect.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def _get_linear_acyl_chain(mol: rdchem.Mol, carbonyl, acyl_neighbor):
    """
    Traverses from the carbonyl carbon through the acyl chain following single bonds between sp3 carbons.
    Returns a list of atoms (starting from the carbonyl) if the chain is strictly linear. Returns None if branching is encountered.
    """
    chain = [carbonyl, acyl_neighbor]
    prev_atom = carbonyl
    current = acyl_neighbor
    # Only follow carbon atoms with sp3 hybridization.
    while True:
        next_atoms = []
        for nbr in current.GetNeighbors():
            if nbr.GetIdx() == prev_atom.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                bond = mol.GetBondBetweenAtoms(current.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    next_atoms.append(nbr)
        if len(next_atoms) == 0:
            break
        elif len(next_atoms) == 1:
            prev_atom = current
            current = next_atoms[0]
            chain.append(current)
        else:
            # Branching encountered
            return None
    return chain

def _atom_heavy_count(mol: rdchem.Mol, atom_idxs):
    """
    Count non-hydrogen atoms among the atoms with indices in atom_idxs.
    """
    count = 0
    for idx in atom_idxs:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 1:
            count += 1
    return count

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule (given as SMILES) is an octanoate ester.
    
    For our purposes an octanoate ester is one where at least one ester group meets all of:
      - The acyl (acid) side is linear and contains exactly eight carbons (the carbonyl carbon plus 7 more),
        with the terminal carbon being a methyl group.
      - The ester oxygen (alcohol side) is not directly or indirectly linked to phosphorus (to filter out, for example, phospholipids).
      - The octanoate fragment (chain + carbonyl oxygen) accounts for a significant fraction of the heavy atoms (>=35%) in the molecule.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): Tuple of (True, reason) if classified as an octanoate ester, else (False, reason).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS to locate an ester: a carbonyl carbon bonded (via a single bond) to an oxygen that is not -OH.
    ester_smarts = "[#6X3](=O)[OX2H0]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return None, None
    
    matches = mol.GetSubstructMatches(ester_pattern)
    if not matches:
        return False, "No ester group found in molecule"
    
    total_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1)
    
    # For each ester candidate, process:
    for match in matches:
        # In the SMARTS, match[0] is the carbonyl carbon, match[1] is the ester oxygen.
        carbonyl = mol.GetAtomWithIdx(match[0])
        ester_oxygen = mol.GetAtomWithIdx(match[1])
        
        # --- Heuristic filter 1: check that the ester oxygen is not part of a phosphate group.
        # Look at the neighbors of the ester oxygen (other than the carbonyl).
        skip_this = False
        for nbr in ester_oxygen.GetNeighbors():
            if nbr.GetIdx() == carbonyl.GetIdx():
                continue
            if nbr.GetAtomicNum() == 15:
                skip_this = True
                break
            # Also check one bond further – if any neighbor of this neighbor is P.
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetAtomicNum() == 15:
                    skip_this = True
                    break
            if skip_this:
                break
        if skip_this:
            continue
        
        # Identify the acyl branch on the carbonyl.
        acyl_neighbors = []
        for nbr in carbonyl.GetNeighbors():
            if nbr.GetIdx() == ester_oxygen.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_neighbors.append(nbr)
        if not acyl_neighbors:
            continue  # no acyl branch found
        
        for acyl in acyl_neighbors:
            # Traverse the acyl chain linearly (starting from the carbonyl)
            chain = _get_linear_acyl_chain(mol, carbonyl, acyl)
            if chain is None:
                continue  # branch encountered
            # For octanoic acid the chain (including the carbonyl carbon) should have exactly 8 carbons.
            if len(chain) != 8:
                continue
            # Check that the terminal carbon is a methyl group (only one carbon neighbor)
            terminal = chain[-1]
            carbon_neighbors = [nbr for nbr in terminal.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                continue
            
            # We consider the octanoate fragment to include the 8-carbon chain plus the carbonyl's double-bonded oxygen.
            # (Find the carbonyl oxygen that is double-bonded.)
            carbonyl_double_oxygen = None
            for nbr in carbonyl.GetNeighbors():
                # Already used ester oxygen; so the other O (if any) should be the double-bonded one.
                if nbr.GetIdx() == ester_oxygen.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
                    if bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        carbonyl_double_oxygen = nbr
                        break
            if carbonyl_double_oxygen is None:
                continue
            
            # Build the set of atom indices for the octanoate fragment.
            fragment_idxs = set([atom.GetIdx() for atom in chain])
            fragment_idxs.add(carbonyl_double_oxygen.GetIdx())
            frag_heavy = _atom_heavy_count(mol, fragment_idxs)
            ratio = frag_heavy / total_heavy if total_heavy > 0 else 0
            
            # We require that the octanoate fragment is a significant fraction (>=35%) of the entire molecule.
            # This heuristic helps filter out cases where an octanoyl chain is just a minor appendage.
            if ratio < 0.35:
                continue
            
            return True, "Molecule contains an octanoate ester moiety (octanoate fragment accounts for {:.0%} of heavy atoms)".format(ratio)
    
    return False, "Molecule does not contain a qualifying octanoate ester moiety"

# Example usage and testing:
if __name__ == "__main__":
    test_smiles_list = [
        # True positives
        "CCCCCCCC(=O)OCC",  # ethyl octanoate
        "CCCCCCCC(=O)OC",   # methyl octanoate
        "CCCCCCCC(=O)OC[C@H](O)CO",  # 1-octanoyl-sn-glycerol
        "CCCCCCCC(=O)OC[C@H](COP(O)(=O)OP(O)(O)=O)OC(=O)CCCCCCCC",  # 1,2-dioctanoyl-sn-glycerol 3-diphosphate
        # Examples that were false positives in previous attempt:
        "P(OC[C@@H](CO)OC(CCCCCCC)=O)(=O)(OCC[N+](C)(C)C)[O-]",  # 2-octanoyl-sn-glycero-3-phosphocholine -> should be False now
        "O1C(C1C/C=C\\CCCCC)CCCCCCCC(O[C@H](COC(=O)CCCCCCC)CO)=O",  # DG(8:0/18:1(12Z)-O(9S,10R)/0:0) -> False
    ]
    for smi in test_smiles_list:
        result, reason = is_octanoate_ester(smi)
        print(f"SMILES: {smi}\nClassification: {result}\nReason: {reason}\n")