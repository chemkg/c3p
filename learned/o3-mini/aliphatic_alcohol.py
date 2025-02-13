"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: Aliphatic Alcohol
Definition: “An alcohol derived from an aliphatic compound.”
An aliphatic alcohol here is taken to be one that has at least one hydroxyl (-OH) group
attached to an sp³-hybridized carbon that is embedded in a predominantly non‐aromatic environment.
This version improves the previous heuristic by:
  • Allowing the –OH-bearing carbon to have sp² (i.e. unsaturated) neighbors if they are not aromatic.
  • If there is only a single heavy neighbor (other than –OH oxygen), and that neighbor is non‐aromatic,
    then we count it as acceptable.
  • Continuing to discount candidates if the carbon is part of a glycoside/sugar ring (i.e. if the ring has >1 oxygen).
"""

from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    
    An aliphatic alcohol is defined as one that has at least one hydroxyl (-OH) group
    attached to an sp³-hybridized carbon. Moreover, the carbon’s local heavy-atom environment 
    (ignoring the hydroxyl oxygen) should be predominantly non‐aromatic.
    Also, if the candidate carbon is part of a ring that looks sugar‐like (has more than one oxygen),
    then that candidate is discounted.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as an aliphatic alcohol, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring info (used to discount sugar-like rings)
    ring_info = mol.GetRingInfo()
    
    candidate_found = False

    # Iterate over all atoms and look for oxygens that may be –OH groups.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # Check if oxygen has at least one hydrogen (detecting –OH groups).
        if atom.GetTotalNumHs() < 1:
            continue
        # We insist that the oxygen is bonded to exactly one heavy (non-H) atom.
        nbrs = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
        if len(nbrs) != 1:
            continue
        neighbor = nbrs[0]
        # We only care if the oxygen is attached to a carbon.
        if neighbor.GetAtomicNum() != 6:
            continue
        # The carbon must be sp³ (this precludes phenolic / vinylic alcohols)
        if neighbor.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue

        # Candidate: we have an –OH on a sp³ carbon.
        # Now inspect the other heavy-neighbors of this carbon (excluding the –OH oxygen).
        heavy_neigh = [n for n in neighbor.GetNeighbors() if n.GetIdx() != atom.GetIdx() and n.GetAtomicNum() > 1]
        # If no other heavy neighbor exists, then (unusual but) allow the candidate.
        if len(heavy_neigh) == 0:
            candidate_found = True
            break

        # Count how many of these neighbors are non-aromatic.
        # Here we count any heavy atom (sp2 or sp3) that is NOT aromatic.
        non_aromatic_count = 0
        total_count = 0
        for nb in heavy_neigh:
            total_count += 1
            if not nb.GetIsAromatic():
                non_aromatic_count += 1

        # If there is only one heavy neighbor and it is non‐aromatic, allow the candidate.
        if len(heavy_neigh) == 1 and non_aromatic_count == 1:
            pass  # candidate qualifies
        elif total_count > 0 and (non_aromatic_count / total_count) < 0.5:
            # If fewer than half of the neighbors are non‐aromatic, skip this candidate.
            continue

        # Check if the candidate carbon is in a ring that looks sugar‐like.
        in_sugar_ring = False
        for ring in ring_info.AtomRings():
            if neighbor.GetIdx() in ring:
                # Count the number of oxygen atoms in this ring.
                oxy_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                if oxy_in_ring > 1:
                    in_sugar_ring = True
                    break
        if in_sugar_ring:
            continue

        # If we have passed all tests, accept this candidate.
        candidate_found = True
        break

    if candidate_found:
        return (True,
                "Molecule contains an aliphatic alcohol group: -OH is attached to a sp³ carbon "
                "in a predominantly non-aromatic (aliphatic) environment.")
    else:
        return (False,
                "No aliphatic alcohol group found based on local environment heuristics.")

# Example usage (for testing purposes)
if __name__ == '__main__':
    test_smiles = [
        "OC(C(C)C)CC/C=C/C=C/C([C@@H](O)C)(C)C",  # Graphostromol G (expected True)
        "CCCCCCCCCCCCC(O)CCCC",                  # nonadecan-5-ol (expected True)
        "CO",                                    # methanol (expected True)
        "O[C@H](C1=CC=C(CO)C=C1)CCCCC",           # 1-[4-(Hydroxymethyl)phenyl]hexan-1-ol (expected True)
        "c1ccccc1O",                            # Phenol (expected False)
        "OCC=C=CC#CC#C",                        # xi-2,3-Octadiene-5,7-diyn-1-ol (previously false negative; expected True)
        "OC(/C=C\\CCCCCCC)C#CC#CC(=O)C=C",       # Falcarinolone (previously false negative; expected True)
        "C(/C=C/CC/C=C\\CC)O",                    # (2E,6Z)-nona-2,6-dien-1-ol (expected True)
    ]
    
    for smi in test_smiles:
        res, reason = is_aliphatic_alcohol(smi)
        print(f"SMILES: {smi}\n  Classified as aliphatic alcohol? {res}\n  Reason: {reason}\n")