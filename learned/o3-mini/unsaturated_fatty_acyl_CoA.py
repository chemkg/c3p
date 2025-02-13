"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: Unsaturated fatty acyl-CoA
Definition: An unsaturated fatty acyl-CoA is formed by the condensation of the thiol group of coenzyme A 
with the carboxyl group of an unsaturated fatty acid. This class is characterized by a thioester linkage 
([C](=O)[S]) connecting an acyl chain to a CoA moiety, where the CoA is recognized by the presence of an adenine 
ring. Furthermore, the acyl chain must be long enough (>= 3 carbons) and contain at least one nonaromatic C=C bond.
"""

from rdkit import Chem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    
    Procedure:
      1. Parse the SMILES into an RDKit molecule; fail if parsing is unsuccessful.
      2. Look for a thioester group ([C](=O)[S]) which links the acyl chain to CoA.
      3. Confirm the presence of a CoA moiety by finding an adenine substructure.
         (We use a more generic pattern 'n1cnc2c1ncnc2' to capture adenine whether or not substituents are present.)
      4. For each thioester match, identify the fatty acyl chain side by finding the neighbor of the carbonyl carbon
         that is not the double-bonded oxygen or sulfur.
      5. Perform a depth-first search (DFS) from that acyl chain anchor to collect connected carbon atoms,
         while avoiding atoms that are part of the CoA (adenine) fragment.
      6. Ensure that the acyl chain contains at least three carbon atoms and at least one nonaromatic C=C bond.
    
    Returns:
      bool:      True if the molecule is classified as an unsaturated fatty acyl-CoA, otherwise False.
      str :      Explanation for the decision.
    """
    # Step 1. Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2. Look for thioester group: [C](=O)[S]
    thioester_pattern = Chem.MolFromSmarts("[C](=O)[S]")
    ts_matches = mol.GetSubstructMatches(thioester_pattern)
    if not ts_matches:
        return False, "Missing thioester group linking the acyl chain to CoA"

    # Step 3. Check for the CoA moiety by identifying the adenine substructure.
    # Updated pattern to match adenine rings even when substitutions are present.
    adenine_pattern = Chem.MolFromSmiles("n1cnc2c1ncnc2")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine substructure; not a recognizable CoA derivative"
    
    adenine_matches = mol.GetSubstructMatches(adenine_pattern)
    adenine_atoms = set()
    for match in adenine_matches:
        adenine_atoms.update(match)

    # Step 4. For each thioester occurrence, try to isolate the fatty acyl chain.
    for match in ts_matches:
        # Each match is a tuple of indices: (carbonyl C, carbonyl O, sulfur)
        carbonyl_idx, oxygen_idx, sulfur_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Identify the acyl-chain attachment (neighbor of the carbonyl which is NOT O or S).
        acyl_anchor = None
        for nb in carbonyl_atom.GetNeighbors():
            if nb.GetIdx() in {oxygen_idx, sulfur_idx}:
                continue
            acyl_anchor = nb
            break
        if acyl_anchor is None:
            continue  # Try next thioester match if this one doesn't yield an acyl chain.
        
        # Step 5. Perform a DFS from acyl_anchor to collect contiguous carbon atoms.
        visited = set()
        stack = [acyl_anchor.GetIdx()]
        chain_atoms = set()
        while stack:
            curr_idx = stack.pop()
            if curr_idx in visited:
                continue
            # Avoid crossing into the CoA moiety (adenine part).
            if curr_idx in adenine_atoms:
                continue
            visited.add(curr_idx)
            atom = mol.GetAtomWithIdx(curr_idx)
            if atom.GetAtomicNum() != 6:  # Only follow carbon atoms.
                continue
            chain_atoms.add(curr_idx)
            # Add neighbors except going back to the carbonyl atom.
            for nb in atom.GetNeighbors():
                if nb.GetIdx() == carbonyl_idx:
                    continue
                if nb.GetIdx() not in visited:
                    stack.append(nb.GetIdx())
        
        # Step 6. Check that the chain is long enough (>= 3 carbons).
        if len(chain_atoms) < 3:
            continue
        
        # Count nonaromatic C=C bonds within the chain.
        unsat_count = 0
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in chain_atoms and a2 in chain_atoms:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and not bond.GetIsAromatic():
                    # Confirm that both atoms are carbon.
                    if (mol.GetAtomWithIdx(a1).GetAtomicNum() == 6 and 
                        mol.GetAtomWithIdx(a2).GetAtomicNum() == 6):
                        unsat_count += 1
        
        # A valid unsaturated fatty acyl-CoA should have at least one nonaromatic double bond.
        if unsat_count >= 1:
            reason = (f"Found a thioester linkage to CoA with an acyl chain of {len(chain_atoms)} carbons "
                      f"and {unsat_count} nonaromatic C=C double bond(s), consistent with an unsaturated fatty acyl-CoA")
            return True, reason

    # If no thioester match produced a qualifying acyl chain:
    return False, "Could not isolate a fatty acyl chain with sufficient length and unsaturation"

# Example usage (when run directly):
if __name__ == "__main__":
    # Testing with one of the provided examples.
    example_smiles = "S(C(=O)CCCCC[C@@H]1[C@@H](C(C=C1)=O)C/C=C\\CC)CCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]2O[C@@H](N3C4=C(C(=NC=N4)N)N=C3)[C@@H]([C@@H]2OP(O)(O)=O)O)(=O)O)(=O)O)(C)C)O)=O"
    result, explanation = is_unsaturated_fatty_acyl_CoA(example_smiles)
    print(result, explanation)