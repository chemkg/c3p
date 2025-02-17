"""
Classifies: CHEBI:25029 leukotriene
"""
"""
Classifies: Leukotrienes – any icosanoid stemming from a C20 polyunsaturated fatty acid derivative generated by leukocytes 
from arachidonic acid. They typically have a 20‐carbon acyclic chain with 4 double bonds (at least 3 conjugated) and a carboxyl group.
The improved heuristic: 1) allow derivatives that are heavier than 500 Da, 2) require the existence of at least one contiguous C20 carbon subchain, 
3) in that subchain require 3 or 4 C=C bonds with ≥3 conjugated, 4) require at least one carboxyl (or related) moiety, and 5) allow at most 1 ring.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_longest_carbon_chain(mol):
    """
    Returns the longest chain (as a list of atom indices) made exclusively of carbon atoms.
    We build an adjacency dictionary on carbon atoms and use DFS to explore paths.
    """
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idxs:
        return []
    neighbors = {idx: [] for idx in carbon_idxs}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i in neighbors and j in neighbors:
            neighbors[i].append(j)
            neighbors[j].append(i)
    longest_chain = []
    def dfs(current, visited, path):
        nonlocal longest_chain
        if len(path) > len(longest_chain):
            longest_chain = path[:]
        for neigh in neighbors[current]:
            if neigh not in visited:
                visited.add(neigh)
                path.append(neigh)
                dfs(neigh, visited, path)
                path.pop()
                visited.remove(neigh)
    for start in carbon_idxs:
        dfs(start, set([start]), [start])
    return longest_chain

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene derivative using improved heuristics:
      1. The molecule must be valid.
      2. It should have at most one ring.
      3. It must have at least one carboxyl group (or related group).
      4. There must exist at least one contiguous subchain (of carbon atoms) of length exactly 20 that
         has 3 or 4 carbon–carbon double bonds (with at least 3 of them conjugated).
      5. (No strict molecular weight filter is applied.)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a leukotriene derivative.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check ring count: leukotrienes are largely acyclic, at most one small ring allowed.
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count > 1:
        return False, f"Too many rings ({ring_count}); leukotrienes are mostly acyclic (at most one small ring allowed)"
    
    # Check for carboxyl group or related moiety.
    # This SMARTS pattern matches carboxylic acids or deprotonated carboxylate groups.
    carboxyl_smarts = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    if not mol.HasSubstructMatch(carboxyl_smarts):
        return False, "No carboxyl (or related) group found; expected in arachidonic acid derivatives"
    
    # Find a long contiguous carbon chain.
    carbon_chain = get_longest_carbon_chain(mol)
    if len(carbon_chain) < 20:
        return False, f"Longest contiguous carbon chain has length {len(carbon_chain)}; a C20 chain is required"
    
    # Check for a subchain of exactly 20 connected carbon atoms meeting the unsaturation criterion.
    # We slide a window of 20 atoms along the longest chain.
    found_valid_subchain = False
    explanation = ""
    for i in range(len(carbon_chain) - 20 + 1):
        subchain = carbon_chain[i:i+20]
        dbl_count = 0
        conjugated_count = 0
        # We look at bonds between consecutive atoms in the chosen subchain.
        for j in range(len(subchain)-1):
            bond = mol.GetBondBetweenAtoms(subchain[j], subchain[j+1])
            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                dbl_count += 1
                if bond.GetIsConjugated():
                    conjugated_count += 1
        # Heuristic: allow either 3 or 4 double bonds in the 20-carbon segment.
        if dbl_count in [3, 4] and conjugated_count >= 3:
            found_valid_subchain = True
            explanation = (f"Found a 20-carbon subchain with {dbl_count} C=C bonds "
                           f"({conjugated_count} of which are conjugated)")
            break
            
    if not found_valid_subchain:
        return False, "No appropriate 20-carbon subchain with the expected unsaturation pattern (3-4 C=C bonds with ≥3 conjugated) was found"
    
    # Build an overall explanation with calculated molecular weight as additional info.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    return True, (f"Molecule has at least one 20-carbon subchain where {explanation}, "
                  f"a carboxyl group is present, MW is {mw:.1f} Da, and ring count is {ring_count} "
                  "– consistent with a leukotriene derivative.")

# Example usage (uncomment to test):
# example_smiles = "CCCCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O"  # leukotriene A4 example
# result, reason = is_leukotriene(example_smiles)
# print(result, reason)