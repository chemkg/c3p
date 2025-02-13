"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: Unsaturated fatty acyl-CoA
Definition: A fatty acyl-CoA is formed from the condensation of the thiol 
group of coenzyme A with the carboxyl group of an unsaturated fatty acid.
This function checks that (1) the molecule can be parsed, (2) it contains a 
thioester linkage characteristic of the acyl-CoA bond, (3) it contains a CoA moiety 
(as evidenced by a specific substructure), and (4) that the “fatty acyl” portion 
(extracted from the thioester linkage after excluding the CoA fragment) is long 
enough and contains at least one nonaromatic C=C bond.
"""

from rdkit import Chem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    
    The approach is:
      1. Parse the SMILES.
      2. Verify a thioester group ([C](=O)[S]) exists.
      3. Verify a CoA moiety exists by matching a characteristic substructure
         (we use "SCCNC(=O)CCNC(=O)" which is common in many CoA molecules).
      4. Locate the acyl chain. From a thioester match, take the carbonyl carbon 
         and then the substituent (neighbor) that is not linked to the sulfur.
      5. Prevent “leakage” into the CoA part by excluding atoms that belong to the CoA match.
      6. Do a simple depth-first search (DFS) starting at that acyl atom
         and “collect” connected aliphatic carbon atoms.
      7. Count the number of carbons (we require at least 5) and count nonaromatic C=C bonds
         among bonds whose both ends are in the acyl fragment.
      
    Returns:
      bool: True if the molecule is judged to be an unsaturated fatty acyl-CoA, False otherwise.
      str : Explanation of the decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the thioester group: look for [C](=O)[S]
    thioester_pattern = Chem.MolFromSmarts("[C](=O)[S]")
    ts_matches = mol.GetSubstructMatches(thioester_pattern)
    if not ts_matches:
        return False, "Missing thioester group linking the fatty acyl part to CoA"
    
    # 2. Check for a CoA moiety.
    # We use a more specific fragment common in CoA: SCCNC(=O)CCNC(=O)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing characteristic CoA moiety"

    # 3. Exclude CoA region from our acyl-chain search.
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    coa_atoms = set()
    for match in coa_matches:
        coa_atoms.update(match)
    
    # 4. For at least one thioester, attempt to extract the fatty acyl chain.
    # For each thioester match, match is a tuple of three atoms: (carbonyl C, carbonyl O, S)
    acyl_chain_atoms = None
    unsat_bonds_count = 0
    chain_carbon_count = 0
    for match in ts_matches:
        carbonyl_idx, oxygen_idx, sulfur_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Determine the acyl-chain attachment: the carbonyl should have two neighbors,
        # one being the oxygen (double bond) and one being the acyl chain (which is not the sulfur).
        neighbors = carbonyl_atom.GetNeighbors()
        acyl_anchor = None
        for nb in neighbors:
            if nb.GetIdx() in {oxygen_idx, sulfur_idx}:
                continue
            # Found the substituent that belongs to the fatty acyl chain.
            acyl_anchor = nb
            break
        if acyl_anchor is None:
            continue  # Try next thioester match if extraction fails.
        
        # 5. Traverse the molecule from the acyl_anchor to collect a contiguous aliphatic chain.
        # We do a DFS from acyl_anchor and collect carbon atoms, but do not enter CoA region.
        visited = set()
        stack = [acyl_anchor.GetIdx()]
        chain_atoms = set()
        while stack:
            idx = stack.pop()
            if idx in visited or idx in coa_atoms:
                continue
            visited.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            # We demand that the atom is carbon.
            if atom.GetAtomicNum() == 6:
                chain_atoms.add(idx)
                # Allow traversal only through non-aromatic C (or if aromatic, we do not count double bonds later)
                for nb in atom.GetNeighbors():
                    # Only continue if neighbor is carbon and not in CoA.
                    if nb.GetAtomicNum() == 6 and nb.GetIdx() not in visited:
                        stack.append(nb.GetIdx())
        # A valid fatty acyl chain tends to have a minimum number of carbons.
        if len(chain_atoms) < 5:
            continue  # Chain too short; try next match.
        
        # 6. Count nonaromatic C=C bonds within the collected chain.
        unsat = 0
        for bond in mol.GetBonds():
            # Check if both atoms belong to the chain and not in CoA region.
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in chain_atoms and a2 in chain_atoms:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and not bond.GetIsAromatic():
                    # Make sure both atoms are carbons.
                    if (mol.GetAtomWithIdx(a1).GetAtomicNum() == 6 and 
                        mol.GetAtomWithIdx(a2).GetAtomicNum() == 6):
                        unsat += 1
        
        # If we found at least one nonaromatic double bond, then we decide this match fulfills our criteria.
        if unsat >= 1:
            acyl_chain_atoms = chain_atoms
            unsat_bonds_count = unsat
            chain_carbon_count = len(chain_atoms)
            break
    
    if acyl_chain_atoms is None:
        return False, "Could not isolate a fatty acyl chain with sufficient length and unsaturation"
    
    # Passed all criteria
    reason = (f"Found a thioester linkage to CoA and an acyl chain with {chain_carbon_count} carbons "
              f"and {unsat_bonds_count} nonaromatic C=C double bond(s), consistent with an unsaturated fatty acyl-CoA")
    return True, reason

# Example usage (if executed as a script):
if __name__ == "__main__":
    # A known positive example (one of the provided true positives)
    example_smiles = "S(C(=O)CCCCC[C@@H]1[C@@H](C(C=C1)=O)C/C=C\\CC)CCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]2O[C@@H](N3C4=C(C(=NC=N4)N)N=C3)[C@@H]([C@@H]2OP(O)(O)=O)O)(=O)O)(=O)O)(C)C)O)=O"
    result, msg = is_unsaturated_fatty_acyl_CoA(example_smiles)
    print(result, msg)