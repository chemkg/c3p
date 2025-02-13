"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: Lipid Hydroperoxide
Defined as any lipid carrying one or more hydroperoxy (-OOH) substituents.
This implementation uses additional constraints:
  1. The molecule must contain at least one -OOH group (detected by the SMARTS "[OX2]-[OX2H]").
  2. The molecule is considered lipid-like if it has:
      • A minimum total number of carbon atoms (>=16).
      • A relatively high molecular weight (>=300 Da).
      • A long uninterrupted carbon chain (at least 8 carbons in a path, ignoring heteroatoms).
      • A limited number of ring systems (at most 1 ring).
If these criteria are met the molecule is classified as a lipid hydroperoxide.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the length of the longest simple path made up exclusively of carbon atoms.
    Builds a graph where nodes are carbon atoms and bonds connect carbons. Then performs
    a depth-first search to find the longest simple (acyclic) chain.
    """
    # Get indices of carbon atoms
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return 0
    # Build a simple graph restricted to carbon atoms.
    carbon_graph = {idx: [] for idx in carbon_atoms}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            carbon_graph[a1.GetIdx()].append(a2.GetIdx())
            carbon_graph[a2.GetIdx()].append(a1.GetIdx())
        
    # DFS to compute longest simple path
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
      3. Have a minimum of 16 carbon atoms.
      4. Possess a molecular weight of at least 300 Da.
      5. Contain a long uninterrupted carbon chain (at least 8 carbons in series).
      6. Have at most one ring system.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a lipid hydroperoxide, else False.
        str: Explanation of the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to detect the hydrogen on the hydroperoxy group.
    mol = Chem.AddHs(mol)
    
    # Identify hydroperoxy group (-OOH) via SMARTS.
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2]-[OX2H]")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy (-OOH) substituent found"
    
    # Count total number of carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 16:
        return False, f"Insufficient number of carbon atoms ({carbon_count}) for a lipid-like structure"
    
    # Check molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, f"Molecular weight too low ({mw:.1f} Da) for a typical lipid"
    
    # Check ring count. (We exclude molecules with more than 1 ring.)
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count > 1:
        return False, f"Too many ring systems ({ring_count}); not typical for a lipid hydroperoxide"
    
    # Determine the length of the longest carbon chain.
    chain_length = longest_carbon_chain(mol)
    if chain_length < 8:
        return False, f"Longest carbon chain is too short ({chain_length} atoms) for a typical lipid tail"
    
    return True, "Molecule contains a hydroperoxy substituent, has a long aliphatic chain, low ring count, and sufficient molecular weight"

# Example usage:
if __name__ == "__main__":
    test_smiles = "CCCC\C=C/C\C=C/C[C@@H](OO)\C=C\C=C/CCCC(O)=O"  # 9-HPETE
    result, reason = is_lipid_hydroperoxide(test_smiles)
    print(result, reason)