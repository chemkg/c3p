"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: Essential Fatty Acid
Definition: Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.
Our heuristic:
  • Find candidate carbonyl anchors from free acid group (SMARTS: "C(=O)[O;H1,-1]") and an acyl ester group (SMARTS: "C(=O)O[C]")
  • For each candidate, follow a linear (non‐branched) walk from the α‐carbon along aliphatic, non‐ring carbons.
    (We count only the alkyl part; the carbonyl carbon is not included in chain length.)
  • Count the alkyl chain’s length and the number of carbon–carbon double bonds.
  • Also compute the fraction of molecule’s total carbon count that appears in that chain.
  • If the alkyl chain is long enough (≥ 12 carbons) and has at least 2 double bonds then:
       – For a free acid candidate, we classify as an essential fatty acid.
       – For an ester candidate, we require that the chain makes up at least 35% of all carbons.
Otherwise the molecule is not classified as an essential fatty acid.
NOTE: This algorithm is heuristic and may still mis‐classify some molecules.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if the molecule (given as a SMILES string) is an essential fatty acid.
    
    The method:
      1. Parse the molecule.
      2. Look for candidate carbonyl groups, from either free acid or ester fragments.
      3. For each candidate, get the α‐carbon(s) (neighbors of the carbonyl carbon that are carbons).
      4. For each α‐carbon candidate, “walk” along a contiguous, non‐branched chain of aliphatic, non‐ring carbons.
         (If there is a branch beyond the first step we stop extension.)
      5. For the extracted chain, count:
            • The alkyl chain length (excluding the carbonyl carbon).
            • The number of C=C bonds between consecutive atoms in the chain.
            • The fraction of all carbons in the full molecule that appear in the chain.
      6. If the chain length is at least 12 and has at least 2 C=C bonds:
            – For free acids (anchor matched by acid SMARTS) we accept immediately.
            – For esters (anchor from ester SMARTS) we accept only if the chain constitutes ≥35% of all carbons.
    
    Args:
      smiles (str): SMILES string for the molecule.
    
    Returns:
      bool: True if the molecule is classified as an essential fatty acid, False otherwise.
      str: Explanation for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Count total number of carbon atoms in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Define SMARTS for a free acid and for a fatty acyl ester.
    acid_smarts = "C(=O)[O;H1,-1]"   # free carboxylic acid group
    ester_smarts = "C(=O)O[C]"         # acyl ester group
    pat_acid = Chem.MolFromSmarts(acid_smarts)
    pat_ester = Chem.MolFromSmarts(ester_smarts)
    
    # Collect candidate anchors.
    # For each match, we record a tuple (anchor_index, candidate_type) where type is 'acid' or 'ester'
    candidate_anchors = []
    for match in mol.GetSubstructMatches(pat_acid):
        # match[0] is the carbonyl carbon.
        candidate_anchors.append((match[0], 'acid'))
    for match in mol.GetSubstructMatches(pat_ester):
        candidate_anchors.append((match[0], 'ester'))
    
    if not candidate_anchors:
        return False, "No carboxyl or acyl ester functional group found."
    
    # Build a dictionary of allowed carbon atoms:
    # Only consider carbon atoms (atomic number 6) that are not in rings.
    allowed = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.IsInRing()):
            idx = atom.GetIdx()
            # List neighbors that are also aliphatic carbons (non-ring).
            nbrs = [nbr.GetIdx() for nbr in atom.GetNeighbors()
                    if nbr.GetAtomicNum() == 6 and not nbr.IsInRing()]
            allowed[idx] = nbrs
    
    best_candidate = None  # to record the best candidate chain that meets criteria
    best_info = None     # tuple: (chain_length, n_double, fraction, candidate_type)
    details = []         # store details about candidate chains for error messaging
    
    # For each candidate anchor, look for adjacent carbons (the α‐carbon(s)).
    for anchor_idx, cand_type in candidate_anchors:
        anchor_atom = mol.GetAtomWithIdx(anchor_idx)
        # For each neighbor of the anchor that is carbon and allowed, consider it as a starting α‐carbon.
        alpha_candidates = [nbr.GetIdx() for nbr in anchor_atom.GetNeighbors() 
                              if nbr.GetAtomicNum() == 6 and (nbr.GetIdx() in allowed)]
        for start in alpha_candidates:
            # Initialize our chain. We will not count the carbonyl carbon (anchor_idx) in the alkyl chain.
            chain = [start]
            current = start
            parent = anchor_idx  # coming from the carbonyl carbon
            # Walk along the chain. For the first step we accept any number of neighbors,
            # but after that we require exactly one allowed extension to enforce linearity.
            while True:
                # Get allowed neighbors of current (only those in allowed) and not equal to parent.
                ext = [n for n in allowed.get(current, []) if n != parent]
                if len(ext) == 1:
                    next_atom = ext[0]
                    chain.append(next_atom)
                    parent, current = current, next_atom
                else:
                    # Either terminal or branching encountered: stop.
                    break

            # Now count chain properties.
            alkyl_length = len(chain)  # chain length, not counting the anchor carbon.
            n_double = 0
            # Count double bonds along the chain (between consecutive atoms in chain).
            for i in range(len(chain)-1):
                bond = mol.GetBondBetweenAtoms(chain[i], chain[i+1])
                if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                    n_double += 1
            # Also check if the terminal atom is a methyl group (preferable for a fatty acid chain)
            term_atom = mol.GetAtomWithIdx(chain[-1])
            # Count non-hydrogen neighbors:
            nonH_neighbors = [nbr for nbr in term_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
            is_methyl = (len(nonH_neighbors) == 1)
            # Compute fraction of molecule carbons in the chain (note: chain does not include the anchor carbon).
            frac = alkyl_length / total_carbons if total_carbons > 0 else 0
            
            candidate_desc = f"Chain from anchor {anchor_idx} (as {cand_type}) with {alkyl_length} carbons, {n_double} double bonds, terminal methyl: {is_methyl}, fraction: {frac:.2f}"
            details.append(candidate_desc)
            
            # Our criteria: need at least 12 alkyl carbons and at least 2 double bonds.
            if alkyl_length >= 12 and n_double >= 2:
                # For ester candidates, require the chain to be a sizeable part of the molecule.
                if cand_type == 'ester' and frac < 0.35:
                    # Candidate chain comes from an ester but the chain is only a small part of a larger molecule.
                    continue
                # Pick candidate with longest chain; if equal, then with more double bonds.
                if best_info is None or (alkyl_length > best_info[0] or (alkyl_length == best_info[0] and n_double > best_info[1])):
                    best_info = (alkyl_length, n_double, frac, cand_type)
                    best_candidate = candidate_desc

    if best_candidate is None:
        return False, "No fatty acyl chain meets criteria: " + " | ".join(details)
    else:
        alkyl_length, n_double, frac, cand_type = best_info
        return True, (f"Contains fatty acyl chain with {alkyl_length} alkyl carbons (excluding carbonyl) and "
                      f"{n_double} C=C bonds (fraction of C atoms: {frac:.2f}) from a {cand_type} group, "
                      "meeting criteria for an essential fatty acid.")

# (Optional) Testing on sample SMILES strings.
if __name__ == '__main__':
    test_smiles = [
        # True positives:
        "CC\\C=C/C\\C=C/C\\C=C/CCCCCC(O)=O",  # all-cis-7,10,13-hexadecatrienoic acid
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCC\\C=C/C\\C=C/C\\C=C/CCCCC",  # phospholipid with qualifying acyl chain
        "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCCCCCCC(O)=O",  # (17Z,20Z,23Z,26Z,29Z)-dotriacontapentaenoic acid
        # False negative (should not be classified):
        "CCCC\\C=C/CCC(O)=O",  # cis-4-decenoic acid → chain too short / too few double bonds
        # False positive example candidates (expected not to pass):
        "OC(=O)CC\\C=C\\C\\C=C\\C\\C=C\\CCCCCCCC",  # 4,7,10,13-Docosatetraenoic acid (here chain extraction might yield a chain that does not dominate)
    ]
    for smi in test_smiles:
        result, reason = is_essential_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")