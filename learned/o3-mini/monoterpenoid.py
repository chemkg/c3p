"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: Monoterpenoid
Definition: Any terpenoid derived from a monoterpene. The term includes compounds in which 
the C10 skeleton of the parent monoterpene has been rearranged or modified by the removal 
of one or more skeletal atoms (generally methyl groups).

Improved Heuristic:
    1. Calculate total carbon and oxygen atoms.
       If the overall oxygen/carbon ratio is high (e.g. >0.4), then the molecule 
       is likely not a typical monoterpenoid.
    2. Extract the Bemis–Murcko scaffold using MurckoScaffold.GetScaffoldForMol.
         a. Count carbons in the scaffold and determine its fraction relative to total carbons.
         b. Check that the scaffold has between 7 and 12 carbon atoms.
         c. Reject the scaffold if more than 50% of its carbons are aromatic.
         d. Also require that the scaffold accounts for at least 50% of the molecule’s total carbon atoms.
    3. If the scaffold is “empty” or does not meet the above criteria, compute the longest 
       continuous carbon chain (using a DFS search) and check that it has between 7 and 12 carbons 
       and represents at least 50% of the molecule’s total carbons.
       
Note: These thresholds (ratio cutoffs, oxygen/carbon ratio etc.) are heuristics that may be tuned.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def longest_carbon_chain(mol):
    """
    Computes the longest continuous, nonrepeating path in the molecule restricted to carbon atoms.
    Returns:
        int: length (number of carbons) of the longest carbon path.
    """
    # Get indices of all carbon atoms
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idxs:
        return 0

    # Build a graph for carbon atoms only.
    graph = {idx: set() for idx in carbon_idxs}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i in graph and j in graph:
            graph[i].add(j)
            graph[j].add(i)

    max_length = 0
    def dfs(node, seen):
        nonlocal max_length
        max_length = max(max_length, len(seen))
        for nb in graph[node]:
            if nb not in seen:
                dfs(nb, seen | {nb})

    for node in graph:
        dfs(node, {node})
    return max_length

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.

    The function uses (a) the Bemis–Murcko scaffold and (b) the longest continuous carbon chain 
    as proxies for the monoterpenoid core. In addition, it checks the fraction of carbons in the 
    core relative to the entire molecule and rejects molecules with an unusually high oxygen/carbon 
    ratio or if the scaffold is predominantly aromatic.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is consistent with a monoterpenoid core, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total carbons and oxygens in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    total_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if total_carbons == 0:
        return False, "No carbon atoms in molecule"
    
    ox_to_c_ratio = total_oxygens / total_carbons
    if ox_to_c_ratio > 0.4:
        return False, f"High oxygen/carbon ratio ({ox_to_c_ratio:.2f}); not consistent with a monoterpenoid"

    # Initialize decision variables.
    accepted = False
    reason = ""
    
    # Attempt to use the Bemis–Murcko scaffold.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is not None:
        scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
        total_scaffold_atoms = scaffold.GetNumAtoms()
        # Count aromatic carbons in scaffold.
        aromatic_scaffold = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
        aromatic_fraction = (aromatic_scaffold / scaffold_carbons) if scaffold_carbons else 0
        
        # Compute fraction of molecule's carbons captured by the scaffold.
        frac_scaffold = scaffold_carbons / total_carbons
        
        if 7 <= scaffold_carbons <= 12:
            if aromatic_fraction > 0.5:
                return False, (f"Scaffold has {scaffold_carbons} carbons but is predominantly aromatic "
                               f"(fraction={aromatic_fraction:.2f}), inconsistent with a monoterpenoid core")
            if frac_scaffold < 0.5:
                return False, (f"Scaffold has {scaffold_carbons} carbons but represents only {frac_scaffold:.2f} "
                               "of total carbons, inconsistent with a monoterpenoid core")
            return True, (f"Scaffold has {scaffold_carbons} carbons, {frac_scaffold:.2f} of molecule, "
                          "consistent with a monoterpenoid core")
        else:
            reason = (f"Scaffold has {scaffold_carbons} carbons, which is not within the expected range (7–12); "
                      f"scaffold represents {frac_scaffold:.2f} of total carbons")
    
    # If scaffold is empty or does not meet criteria, check the longest continuous carbon chain.
    longest_chain = longest_carbon_chain(mol)
    frac_chain = longest_chain / total_carbons
    if 7 <= longest_chain <= 12:
        if frac_chain < 0.5:
            return False, (f"Longest carbon chain has {longest_chain} carbons but represents only "
                           f"{frac_chain:.2f} of total carbons, inconsistent with a monoterpenoid core")
        return True, (f"Longest carbon chain has {longest_chain} carbons, representing {frac_chain:.2f} of molecule, "
                      "consistent with a monoterpenoid core")
    else:
        if reason:
            reason += f"; additionally, longest carbon chain has {longest_chain} carbons"
        else:
            reason = f"Longest carbon chain has {longest_chain} carbons, which is not within the expected range (7–12)"
    
    return False, reason

# Example usage (for testing):
if __name__ == "__main__":
    test_examples = [
        # True positives
        "O(C(CCC=C(C)C)(C)C=C)C(=O)CCCCCCC",  # Linalyl caprylate
        "CC1(C)C2CCC1(C=O)C(=O)C2",            # Oxocamphor
        "OCCCC1(C2(C3C2CC1C3)C)",               # Tricycloekasantalol
        "OC(=O)[C@@H]1[C@]2(O)[C@]([C@H](OC(=O)C)C([C@@H]2[C@@H](O)C[C@H]1C)(C)C)(CO)C",  # methyl 7alpha-acetoxydeacetylbotryoloate
        "C(CC=C(C)C)/C(/C)=C/C",               # 2,6-Octadiene, 2,6-dimethyl-
        # False positive example from previous attempt:
        "[H]C(=CC(=O)C(O)=O)C(CC(O)=O)C(O)=O"  # 5-oxopent-3-ene-1,2,5-tricarboxylic acid
    ]
    for sm in test_examples:
        res, msg = is_monoterpenoid(sm)
        print(f"SMILES: {sm}\nResult: {res}\nReason: {msg}\n")