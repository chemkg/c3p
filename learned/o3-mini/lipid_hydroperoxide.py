"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: Lipid Hydroperoxide
Defined as any lipid carrying one or more hydroperoxy (-OOH) substituents.
This implementation uses additional constraints:
  1. The molecule must contain at least one hydroperoxy (-OOH) group (detected by SMARTS "[OX2H]-[OX2]").
  2. The molecule must be lipid‐like:
       • At least 16 carbon atoms.
       • Molecular weight of at least 300 Da.
       • A long uninterrupted carbon chain of at least 8 carbons.
       • The longest carbon chain must account for at least 50% of all carbons.
       • At most one ring system.
       • An oxygen-to-carbon ratio below 0.3.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the length of the longest simple path made exclusively of carbon atoms.
    Builds a graph where only carbon atoms (atomic number 6) are nodes and bonds connecting them are edges.
    A DFS is used to find the longest simple (acyclic) chain.
    """
    # Identify carbon atoms
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return 0
    # Build a graph (dictionary) restricted to carbon atoms.
    carbon_graph = {idx: [] for idx in carbon_atoms}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            carbon_graph[a1.GetIdx()].append(a2.GetIdx())
            carbon_graph[a2.GetIdx()].append(a1.GetIdx())
    
    # Depth-first search to find longest simple path
    def dfs(node, visited):
        max_length = 1
        for neigh in carbon_graph[node]:
            if neigh not in visited:
                length = 1 + dfs(neigh, visited | {neigh})
                if length > max_length:
                    max_length = length
        return max_length
    
    max_chain = 0
    for node in carbon_graph:
        chain_length = dfs(node, {node})
        if chain_length > max_chain:
            max_chain = chain_length
    return max_chain

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    
    The molecule must:
      1. Parse correctly.
      2. Contain at least one hydroperoxy (-OOH) group.
      3. Have at least 16 carbon atoms.
      4. Possess a molecular weight of at least 300 Da.
      5. Have a long uninterrupted carbon chain (>=8 atoms) that also accounts for at least 50% of all carbons.
      6. Have an oxygen:carbon ratio of less than 0.3.
      7. Have at most one ring system.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a lipid hydroperoxide, else False.
      str: Explanation of the classification.
    """
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (to capture hydroperoxy H properly)
    mol = Chem.AddHs(mol)
    
    # Check for the hydroperoxy group (-OOH) using SMARTS:
    # The pattern "[OX2H]-[OX2]" matches an -OOH group.
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2H]-[OX2]")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy (-OOH) substituent found"
    
    # Count carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 16:
        return False, f"Insufficient carbon count: {carbon_count}; need at least 16 for lipid-like structure"
    
    # Check molecular weight
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, f"Molecular weight too low ({mw:.1f} Da) for a typical lipid"
    
    # Check ring systems (allow at most one)
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count > 1:
        return False, f"Too many ring systems ({ring_count}); not typical for a lipid hydroperoxide"
    
    # Determine the longest uninterrupted carbon chain length
    chain_length = longest_carbon_chain(mol)
    if chain_length < 8:
        return False, f"Longest carbon chain too short ({chain_length} atoms); need at least 8"
    
    # Check that the main chain represents a large fraction of the molecule’s carbons
    if chain_length < 0.5 * carbon_count:
        return False, f"Longest carbon chain ({chain_length} atoms) represents less than 50% of total carbons ({carbon_count})"
    
    # Compute oxygen count and check oxygen:carbon ratio
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    oxygen_ratio = oxygen_count / carbon_count if carbon_count > 0 else 0
    if oxygen_ratio >= 0.3:
        return False, f"Oxygen to carbon ratio too high ({oxygen_ratio:.2f}); not typical for a lipid hydroperoxide"
    
    return True, "Molecule contains a hydroperoxy substituent, a long main carbon chain, low ring count, and acceptable oxygen content for a lipid hydroperoxide"

# Example usage:
if __name__ == "__main__":
    # Test with one known lipid hydroperoxide: 9-HPETE
    test_smiles = "CCCC\\C=C/C\\C=C/CC(OO)\\C=C\\C=C/CCCC(O)=O"
    result, reason = is_lipid_hydroperoxide(test_smiles)
    print(result, reason)