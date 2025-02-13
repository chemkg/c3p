"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: Lipid Hydroperoxide
Defined as any lipid carrying one or more hydroperoxy (-OOH) substituents.
This revised implementation now checks:
  1. The molecule parses correctly.
  2. There is at least one hydroperoxy (-OOH) group (using SMARTS "[OX2H]-[OX2]").
  3. The molecule has at least 16 carbon atoms and a molecular weight of at least 300 Da.
  4. It has an uninterrupted carbon chain of at least 8 atoms that represents at least 50% of total carbons.
  5. It has at most one ring system.
  6. Its oxygen-to-carbon ratio is below 0.35.
  7. It does not contain anionic oxygen atoms that are not part of a carboxylate group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the length of the longest simple (acyclic) carbon chain in the molecule.
    Only carbon atoms (atomic number 6) and bonds among them are considered.
    Uses a DFS approach and returns the maximum chain length found.
    """
    # Get indexes of carbon atoms
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return 0
    # Build a graph (dictionary) for carbon atoms only.
    carbon_graph = {idx: [] for idx in carbon_atoms}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            carbon_graph[a1.GetIdx()].append(a2.GetIdx())
            carbon_graph[a2.GetIdx()].append(a1.GetIdx())
    
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
      5. Have a long uninterrupted carbon chain (>=8 atoms) that accounts for at least 50% of all carbons.
      6. Have an oxygen:carbon ratio below 0.35.
      7. Have at most one ring system.
      8. Not contain any anionic oxygen atoms that are not part of a carboxylate group.
      
    Args:
      smiles (str): SMILES string of the molecule.

    Returns:
      bool: True if the molecule is a lipid hydroperoxide; False otherwise.
      str: Explanation of the decision.
    """
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydroperoxy H is explicit
    mol = Chem.AddHs(mol)
    
    # Check for hydroperoxy group (-OOH) using SMARTS pattern.
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2H]-[OX2]")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy (-OOH) substituent found"
    
    # Count carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 16:
        return False, f"Too few carbon atoms ({carbon_count}); require at least 16 for a lipid-like structure"
    
    # Check molecular weight
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, f"Molecular weight too low ({mw:.1f} Da); require at least 300 Da for a typical lipid"
    
    # Check ring systems (allow at most one)
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count > 1:
        return False, f"Too many ring systems ({ring_count}); not typical for a lipid hydroperoxide"
    
    # Check for a long uninterrupted carbon chain
    chain_length = longest_carbon_chain(mol)
    if chain_length < 8:
        return False, f"Longest carbon chain is too short ({chain_length} atoms); need at least 8"
    if chain_length < 0.5 * carbon_count:
        return False, f"Longest carbon chain ({chain_length} atoms) is less than 50% of total carbons ({carbon_count})"
    
    # Compute oxygen count and check oxygen-to-carbon ratio (< 0.35)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    oxygen_ratio = oxygen_count / carbon_count if carbon_count > 0 else 0
    if oxygen_ratio >= 0.35:
        return False, f"Oxygen to carbon ratio too high ({oxygen_ratio:.2f}); not typical for a lipid hydroperoxide"
    
    # Exclude molecules with anionic oxygens that are not part of a carboxylate group.
    # First, collect all oxygen atoms with formal charge -1.
    anionic_oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1]
    
    # Define a SMARTS pattern for a carboxylate group: [CX3](=O)[O-]
    carboxylate_smarts = Chem.MolFromSmarts("[$([CX3](=O)[O-])]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_smarts)
    carboxylate_oxygens = set()
    for match in carboxylate_matches:
        # Add indices of oxygen atoms in the match that have a -1 charge.
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
                carboxylate_oxygens.add(idx)
    
    # If any anionic oxygen is present that is not part of a carboxylate, then reject.
    for atom in anionic_oxygens:
        if atom.GetIdx() not in carboxylate_oxygens:
            return False, "Molecule contains anionic oxygen groups not typical for a lipid hydroperoxide"
    
    return True, "Molecule contains a hydroperoxy substituent, a long main carbon chain, low ring count, and acceptable oxygen content for a lipid hydroperoxide"

# Example usage:
if __name__ == "__main__":
    # Test with one known lipid hydroperoxide: 9-HPETE
    test_smiles = "CCCC\\C=C/C\\C=C/CC(OO)\\C=C\\C=C/CCCC(O)=O"
    result, reason = is_lipid_hydroperoxide(test_smiles)
    print(result, reason)