"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: Aliphatic Alcohol 
Definition: “An alcohol derived from an aliphatic compound”
An aliphatic alcohol here is taken to be one that contains at least one hydroxyl (-OH) group 
attached to an sp3-hybridized carbon that is embedded in a predominantly aliphatic (saturated) environment.
We use a two‐step procedure:
  1. Identify –OH groups attached to sp³ carbons.
  2. For each candidate, look at the local atomic neighborhood of that carbon (ignoring the –OH oxygen) 
     and compute the fraction of neighbors that are both sp³ and non‐aromatic.
     Also, if the candidate carbon is part of a ring that appears to be part of a glycoside/sugar (i.e. has more than one oxygen),
     we discount that candidate.
If at least one candidate passes these tests then the molecule is classified as an aliphatic alcohol.
Otherwise it is not.
"""

from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    
    An aliphatic alcohol is defined here as one that has at least one hydroxyl (-OH) group 
    attached to an sp³-hybridized carbon that is in a predominantly aliphatic (non‐aromatic, saturated)
    environment. In addition, we discount cases (such as glycosidic rings) where the –OH is part of a cyclic sugar‐like 
    system.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an aliphatic alcohol, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ri = mol.GetRingInfo()  # ring information which we will use to flag sugar-like rings

    candidate_found = False

    # Iterate over all atoms looking for oxygen atoms that might be –OH groups.
    for atom in mol.GetAtoms():
        # Look for oxygen (atomic number 8)
        if atom.GetAtomicNum() != 8:
            continue
        # Check if this oxygen is bonded to at least one hydrogen.
        # (GetTotalNumHs() counts both implicit and explicit hydrogens.)
        if atom.GetTotalNumHs() < 1:
            continue
        
        # We further require that it is attached to exactly one heavy neighbor.
        # (This avoids confusing –OH in esters or carboxyl groups with genuine alcohol groups.)
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 1:
            continue
            
        neighbor = neighbors[0]
        # We want the –OH to hang off a carbon atom.
        if neighbor.GetAtomicNum() != 6:
            continue
        
        # The candidate carbon must be sp3 (this avoids benzylic or vinylic alcohols)
        if neighbor.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue

        # At this point we have a candidate alcohol group.
        # Rather than immediately accepting it, we check its local environment.
        # We look at the neighbors of the alcohol-bearing carbon excluding the –OH oxygen.
        nbrs = [a for a in neighbor.GetNeighbors() if a.GetIdx() != atom.GetIdx()]
        if len(nbrs) == 0:
            # Unusual case, but if the carbon has no other neighbor then let it pass.
            candidate_found = True
            break

        # Calculate the fraction of neighbors that are both sp3-hybridized and non-aromatic.
        total = 0
        count_aliphatic = 0
        for nb in nbrs:
            # Consider only heavy atoms.
            if nb.GetAtomicNum() == 1:
                continue
            total += 1
            if (nb.GetHybridization() == Chem.rdchem.HybridizationType.SP3) and (not nb.GetIsAromatic()):
                count_aliphatic += 1
        # If the fraction is too low, then the alcohol-bearing carbon is in a mixed or non-aliphatic environment.
        if total > 0 and (count_aliphatic / total) < 0.5:
            continue

        # Next, if the candidate carbon is part of a ring, check whether the ring looks sugar-like,
        # e.g. if the ring has more than one oxygen then we discount it.
        in_ring = False
        for ring in ri.AtomRings():
            if neighbor.GetIdx() in ring:
                # Count oxygens in this ring
                oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                # For a simple aliphatic ring, we expect few oxygens.
                # If more than one oxygen is present, the ring likely belongs to a glycoside or similar.
                if oxy_count > 1:
                    in_ring = True
                    break
        if in_ring:
            continue

        # If we get here, we have a valid candidate.
        candidate_found = True
        break

    if candidate_found:
        return True, "Molecule contains an aliphatic alcohol group: -OH is attached to an sp3 carbon in an aliphatic environment."
    else:
        return False, "No aliphatic alcohol group found based on local environment heuristics."

# Example usage (you can comment out or remove these when using as a module)
if __name__ == '__main__':
    # A few simple examples including one known false-negative test case:
    examples = [
        "OC(C(C)C)CC/C=C/C=C/C([C@@H](O)C)(C)C",  # Graphostromol G (expected True)
        "CCCCCCCCCCCCC(O)CCCC",                  # nonadecan-5-ol (expected True)
        "CO",                                    # methanol (expected True)
        "O[C@H](C1=CC=C(CO)C=C1)CCCCC",           # 1-[4-(Hydroxymethyl)phenyl]hexan-1-ol (expected True, even though -OH carbon is attached to aromatic ring on one side)
        "c1ccccc1O"                             # Phenol (expected False – OH on aromatic carbon)
    ]
    
    for smi in examples:
        res, reason = is_aliphatic_alcohol(smi)
        print(f"SMILES: {smi}\n  Classified as aliphatic alcohol? {res}\n  Reason: {reason}\n")