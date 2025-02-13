"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: Tetraterpenoid 
Defined as any terpenoid derived from a tetraterpene (typically with a C40 core,
possibly rearranged or slightly modified), which commonly exhibits an extended
conjugated polyene chain. This improved version looks for a long contiguous conjugated
chain (of carbon atoms) in the molecule.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Improved classification consists of:
      1. Finding the longest connected chain of carbon atoms where the connecting bonds are conjugated.
         This is taken as a proxy for the tetraterpene-derived conjugated backbone.
      2. The longest backbone must roughly contain between 30 and 50 carbon atoms.
      3. Along that backbone, the number of non‐aromatic C=C bonds should be at least 7,
         consistent with an extended conjugated polyene system.
      4. The overall molecular weight must be within the range 300–750 Da.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a tetraterpenoid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # First, check overall molecular weight
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300 or mw > 750:
        return False, f"Molecular weight ({mw:.1f}) is not in the expected range (300-750 Da) for tetraterpenoids."
    
    # Build a carbon-backbone graph using only bonds that are conjugated.
    # Only include carbon atoms (atomic number 6) that participate in at least one conjugated bond.
    backbone_graph = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            idx = atom.GetIdx()
            backbone_graph[idx] = []  # initialize neighbor list

    for bond in mol.GetBonds():
        # We take any bond that is flagged as conjugated.
        if bond.GetIsConjugated():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Only include if both atoms are carbons
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                i1 = a1.GetIdx()
                i2 = a2.GetIdx()
                # Only add if the atoms are in our backbone_graph
                if i1 in backbone_graph and i2 in backbone_graph:
                    backbone_graph[i1].append(i2)
                    backbone_graph[i2].append(i1)
    
    # If we found no conjugated carbon bonds, this is not a tetraterpenoid.
    if not backbone_graph:
        return False, "No carbon atoms available for conjugated backbone search."
    
    # Now, we try to find the longest simple path in the backbone graph.
    # Because the graphs are small (dozens of atoms), a brute-force DFS is acceptable.
    def dfs(current, visited):
        best_path = [current]
        for neighbor in backbone_graph[current]:
            if neighbor not in visited:
                candidate = [current] + dfs(neighbor, visited | {neighbor})
                if len(candidate) > len(best_path):
                    best_path = candidate
        return best_path

    longest_path = []
    for node in backbone_graph:
        path = dfs(node, {node})
        if len(path) > len(longest_path):
            longest_path = path

    backbone_length = len(longest_path)
    # Check if the longest conjugated carbon backbone is in the expected range.
    if backbone_length < 30 or backbone_length > 50:
        return False, (f"Longest conjugated carbon chain length ({backbone_length}) is not "
                       "in the expected range (30-50) for tetraterpenoids.")
    
    # Now, count the number of non-aromatic double bonds along this backbone.
    # For each pair of consecutive atoms in the longest path, count the bond if it is a non-aromatic double bond.
    double_bond_count = 0
    for i in range(len(longest_path) - 1):
        bond = mol.GetBondBetweenAtoms(longest_path[i], longest_path[i+1])
        if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
            if not (bond.GetBeginAtom().GetIsAromatic() or bond.GetEndAtom().GetIsAromatic()):
                double_bond_count += 1

    if double_bond_count < 7:
        return False, (f"Found only {double_bond_count} non-aromatic C=C bonds along the backbone; "
                       "expected at least 7 for an extended conjugated system.")

    return True, (f"Longest conjugated carbon chain has {backbone_length} carbons with "
                  f"{double_bond_count} non-aromatic C=C bonds, and molecular weight ({mw:.1f}) "
                  "is in the expected range for a tetraterpenoid.")

# Example usage:
# Uncomment the lines below to test with one of the example SMILES.
# test_smiles = "OC1CC(C(=C(C1=O)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=C(C)C)\\C)C"
# result, reason = is_tetraterpenoid(test_smiles)
# print(result, reason)