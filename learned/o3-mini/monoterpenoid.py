"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: Monoterpenoid
Definition: Any terpenoid derived from a monoterpene. The term includes compounds in which 
the C10 skeleton of the parent monoterpene has been rearranged or modified by the removal 
of one or more skeletal atoms (generally methyl groups).

Heuristic:
    1. Extract the Bemis–Murcko scaffold using MurckoScaffold.GetScaffoldForMol.
       If the scaffold is nonempty, count the number of carbon atoms in it.
         - Also, if the scaffold is mainly aromatic (for example, diphenyl ethers) then reject.
    2. If the scaffold is “empty” (often for acyclic terpenoids) then compute the longest continuous 
       chain of carbon atoms in the molecule and use that number.
    3. Accept if the number (scaffold carons or longest chain) is in the range 7–12.
    
Note: This heuristic still may mis-classify ambiguous cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def longest_carbon_chain(mol):
    """
    Computes the longest continuous path (without repeating atoms) 
    in the molecule restricted to carbon atoms.
    """
    # Get all carbon atom indices
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idxs:
        return 0
    
    # Build a graph (dictionary) of carbon atoms: key is atom index, value is a set of neighbor indices (only carbons).
    graph = {idx: set() for idx in carbon_idxs}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i in graph and j in graph:  # both carbon atoms
            graph[i].add(j)
            graph[j].add(i)
    
    # DFS to find the longest simple path in this undirected graph.
    max_length = 0
    # We'll use recursion with a seen set; since molecules are small the exponential cost is acceptable.
    def dfs(node, seen):
        nonlocal max_length
        max_length = max(max_length, len(seen))
        for nb in graph[node]:
            if nb not in seen:
                dfs(nb, seen | {nb})
    
    # Try starting from every carbon
    for node in graph:
        dfs(node, {node})
    
    return max_length

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.

    The function first attempts to extract the molecule's Bemis–Murcko scaffold (core framework)
    and counts its carbon atoms. If the scaffold has a carbon count within the range of 7–12 and is not 
    heavily aromatic, it is considered consistent with a monoterpenoid core.
    
    For acyclic cases where the Murcko scaffold is empty (or has zero carbons), the function computes 
    the longest continuous chain of carbons in the molecule and checks if that chain length is between 7 and 12.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is consistent with a monoterpenoid core, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Extract the Bemis–Murcko scaffold from the molecule
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        scaffold_carbons = 0
    else:
        scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Use the scaffold if it has carbons; if not, use the longest carbon chain metric.
    if scaffold_carbons > 0:
        # Check if the scaffold’s carbon count is in the accepted range.
        if 7 <= scaffold_carbons <= 12:
            # Also, check the fraction of aromatic carbons.
            aromatic_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
            aromatic_fraction = aromatic_carbons / scaffold_carbons if scaffold_carbons > 0 else 0
            if aromatic_fraction > 0.5:
                return False, f"Scaffold has {scaffold_carbons} carbons but is predominantly aromatic (fraction={aromatic_fraction:.2f}), inconsistent with a monoterpenoid core"
            else:
                return True, f"Scaffold has {scaffold_carbons} carbons, consistent with a monoterpenoid core"
        else:
            return False, f"Scaffold has {scaffold_carbons} carbons, which is not consistent with a monoterpenoid core (expected 7–12 carbons)"
    else:
        # For acyclic molecules the Murcko scaffold may be empty.
        longest_chain = longest_carbon_chain(mol)
        if 7 <= longest_chain <= 12:
            return True, f"Longest carbon chain has {longest_chain} carbons, consistent with a monoterpenoid core"
        else:
            return False, f"Longest carbon chain has {longest_chain} carbons, which is not consistent with a monoterpenoid core (expected 7–12 carbons)"

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        "O(C(CCC=C(C)C)(C)C=C)C(=O)CCCCCCC",  # Linalyl caprylate (acyclic monoterpenoid ester)
        "CC1(C)C2CCC1(C=O)C(=O)C2",            # Oxocamphor (cyclic)
        "OCCCC1(C2(C3C2CC1C3)C)",               # Tricycloekasantalol
        "S1N(C2=CC=C(C=C2)C(OCC)=O)C(=O)C=C1"   # A false positive in previous attempt
    ]
    for sm in test_smiles:
        result, reason = is_monoterpenoid(sm)
        print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")