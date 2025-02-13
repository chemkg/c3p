"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: Unsaturated fatty acyl-CoA
Definition: An unsaturated fatty acyl-CoA is produced from the condensation of the
thiol group of coenzyme A with the carboxyl group of an unsaturated fatty acid.
This function checks that (1) the molecule can be parsed; (2) it contains a thioester
group linking an acyl chain with the CoA moiety; (3) it contains evidence of CoA by
finding an adenine ring; and (4) from the thioester we extract the fatty acyl chain and
verify it is long enough (>=3 carbons) and contains at least one nonaromatic C=C bond.
"""

from rdkit import Chem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.

    The procedure is as follows:
      1. Parse the SMILES into an RDKit molecule. Fail if the SMILES is invalid.
      2. Look for a thioester group (SMARTS: [C](=O)[S]). This is required for acyl-CoA.
      3. Check that a CoA moiety is present by finding an adenine substructure
         (SMARTS: n1cnc2c(n1)nc(n2)) which is present in all coenzyme A derivatives.
      4. For each thioester match, identify the fatty acyl side. The carbonyl carbon
         should have two neighbors, one being the double-bonded oxygen and the other
         being the acyl chain (the one not attached to the sulfur).
      5. From that acyl-chain “anchor,” perform a depth-first search (DFS) to collect
         contiguous carbon atoms. To avoid “leakage” into the CoA part we skip any
         atom that is part of the adenine match.
      6. Count the number of carbons in the acyl chain and the number of nonaromatic
         C=C double bonds (only those bonds whose both ends lie in the collected chain).
      7. If the chain has at least three carbon atoms and at least one nonaromatic C=C bond,
         declare the molecule to belong to the unsaturated fatty acyl-CoA class.

    Returns:
      bool: True if the molecule is considered an unsaturated fatty acyl-CoA, otherwise False.
      str : Explanation of the decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # 1. Check for thioester group [C](=O)[S]
    thioester_pattern = Chem.MolFromSmarts("[C](=O)[S]")
    ts_matches = mol.GetSubstructMatches(thioester_pattern)
    if not ts_matches:
        return False, "Missing thioester group linking the acyl chain to CoA"
        
    # 2. Check for a CoA motif by looking for the adenine substructure.
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(n1)nc(n2)")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine substructure; not a recognizable CoA derivative"

    adenine_matches = mol.GetSubstructMatches(adenine_pattern)
    adenine_atoms = set()
    for match in adenine_matches:
        adenine_atoms.update(match)
        
    # 3. For each thioester group, try to isolate the fatty acyl chain.
    # We loop over thioester matches and try to extract the acyl chain.
    for match in ts_matches:
        # match is a tuple of indices: (carbonyl C, carbonyl O, sulfur)
        carbonyl_idx, oxygen_idx, sulfur_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Identify the acyl-chain attachment: the carbonyl should have two neighbors:
        # one is the double-bonded oxygen and one is the acyl chain (the one not being S)
        acyl_anchor = None
        for nb in carbonyl_atom.GetNeighbors():
            if nb.GetIdx() in {oxygen_idx, sulfur_idx}:
                continue
            # This neighbor is taken as the start of the acyl chain.
            acyl_anchor = nb
            break
        # If no suitable acyl chain found, try the next thioester match.
        if acyl_anchor is None:
            continue

        # 4. From the acyl_anchor, perform a DFS to collect contiguous carbon atoms.
        # We only follow atoms that are carbon and we do not allow entering atoms that
        # are part of the adenine (CoA) substructure.
        visited = set()
        stack = [acyl_anchor.GetIdx()]
        chain_atoms = set()
        while stack:
            curr_idx = stack.pop()
            if curr_idx in visited:
                continue
            # Do not cross into atoms that are clearly part of the CoA (adenine)
            if curr_idx in adenine_atoms:
                continue
            visited.add(curr_idx)
            atom = mol.GetAtomWithIdx(curr_idx)
            if atom.GetAtomicNum() != 6:
                continue  # Only follow carbon atoms
            chain_atoms.add(curr_idx)
            # Add neighbors (except going back to the carbonyl atom)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() == carbonyl_idx:
                    continue
                if nb.GetIdx() not in visited:
                    stack.append(nb.GetIdx())
            
        # 5. Check that the chain is long enough.
        if len(chain_atoms) < 3:
            continue  # Too short: try next thioester match.
            
        # 6. Count nonaromatic C=C bonds within the collected chain.
        unsat_count = 0
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in chain_atoms and a2 in chain_atoms:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and not bond.GetIsAromatic():
                    # Both atoms should be carbons.
                    if (mol.GetAtomWithIdx(a1).GetAtomicNum() == 6 and 
                        mol.GetAtomWithIdx(a2).GetAtomicNum() == 6):
                        unsat_count += 1
        # If at least one nonaromatic double bond is found, classify as unsaturated.
        if unsat_count >= 1:
            reason = (f"Found a thioester linkage to CoA with an acyl chain of {len(chain_atoms)} carbons "
                      f"and {unsat_count} nonaromatic C=C double bond(s), consistent with an unsaturated fatty "
                      "acyl-CoA")
            return True, reason

    # If no thioester match yielded a qualifying fatty acyl chain:
    return False, "Could not isolate a fatty acyl chain with sufficient length and unsaturation"

# Example usage (when run directly):
if __name__ == "__main__":
    # Try one of the provided positive examples:
    example_smiles = "S(C(=O)CCCCC[C@@H]1[C@@H](C(C=C1)=O)C/C=C\\CC)CCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]2O[C@@H](N3C4=C(C(=NC=N4)N)N=C3)[C@@H]([C@@H]2OP(O)(O)=O)O)(=O)O)(=O)O)(C)C)O)=O"
    result, explanation = is_unsaturated_fatty_acyl_CoA(example_smiles)
    print(result, explanation)