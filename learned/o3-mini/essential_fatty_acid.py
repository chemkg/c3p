"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: Essential Fatty Acid
Definition: Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.
This improved function tries all candidate carbonyl carbons (from a free acid or ester fatty acyl group),
and for each candidate it looks for a contiguous, linear chain of non‐ring carbons.
After extracting the chain, it counts the number of carbon atoms (chain length) and the number of carbon–carbon double bonds.
If the chain is long enough (≥ 12 carbons) and has at least 2 C=C bonds, the molecule will be classified as an essential fatty acid.

Note: This heuristic “chain‐extraction” approach may still mis‐classify some molecules.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    The method is as follows:
      1. Verify the SMILES string is valid.
      2. Look for candidate functional groups: a free carboxylic acid (SMARTS: "C(=O)[O;H1,-1]")
         or a fatty acyl ester (SMARTS: "C(=O)O[C]").
      3. For each candidate carbonyl carbon (the “anchor”), find its neighboring carbons that could be
         part of a fatty acyl chain.
      4. For each starting neighbor, “walk” along a contiguous chain of aliphatic, non‐ring carbons.
         After the first extension (which is allowed to be ambiguous), every subsequent step is only made
         if there is exactly one way to continue (a strictly “linear” extension).
      5. Count the total number of carbons in the chain and tally the number of double bonds among consecutive carbons.
      6. If for any candidate chain the chain length is at least 12 and the number of C=C bonds is at least 2,
         then classify the molecule as an essential fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an essential fatty acid, else False.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define SMARTS for a free acid and for an ester group.
    acid_smarts = "C(=O)[O;H1,-1]"   # free carboxylic acid group 
    ester_smarts = "C(=O)O[C]"         # acyl ester group (e.g., in phospholipids)
    pat_acid = Chem.MolFromSmarts(acid_smarts)
    pat_ester = Chem.MolFromSmarts(ester_smarts)
    
    # Collect candidate anchor indices (the carbonyl carbon in the acid/ester).
    candidate_indices = set()
    for match in mol.GetSubstructMatches(pat_acid):
        # match[0] is the carbonyl carbon.
        candidate_indices.add(match[0])
    for match in mol.GetSubstructMatches(pat_ester):
        candidate_indices.add(match[0])
        
    if not candidate_indices:
        return False, "No carboxylic acid (or fatty acyl ester) functional group found."
    
    # Build a dictionary of "allowed" carbon atoms.
    # Only include carbon atoms (atomic number 6) that are not in a ring.
    # These atoms are candidates for being part of a linear acyl chain.
    allowed = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            idx = atom.GetIdx()
            # List neighbors (only carbons, not in rings) by index.
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and not nbr.IsInRing()]
            allowed[idx] = neighbors

    # To record the best fatty acyl chain found.
    best_chain = None
    best_chain_info = (0, 0)  # (chain length, number of double bonds)
    
    # For each candidate anchor, attempt to extract a fatty acyl chain.
    # The anchor is the carbonyl carbon. We look at its neighbors that are carbons and in our allowed graph.
    for anchor in candidate_indices:
        # Get neighbor carbons (skip any non-carbon neighbors)
        anchor_atom = mol.GetAtomWithIdx(anchor)
        neighbor_candidates = [nbr.GetIdx() for nbr in anchor_atom.GetNeighbors() 
                               if nbr.GetAtomicNum() == 6 and (nbr.GetIdx() in allowed)]
        # Try for each neighbor as starting point of a chain.
        for start in neighbor_candidates:
            # Initialize the chain with the anchor and the starting neighbor.
            chain = [anchor, start]
            # To extend the chain, use a simple linear “walk.”
            # After the first step (which may have multiple options) we require that the extension is unique.
            current = start
            parent = anchor
            while True:
                # Get allowed neighbors of the current carbon (excluding the one we just came from).
                ext = [n for n in allowed.get(current, []) if n != parent]
                # If exactly one way to extend, continue the chain.
                if len(ext) == 1:
                    next_atom = ext[0]
                    chain.append(next_atom)
                    parent, current = current, next_atom
                else:
                    break  # either a terminal atom or a branch point, so stop extending.
            # Count chain length (number of carbons in our extracted chain).
            chain_length = len(chain)
            # Also, count the number of C=C double bonds along the chain.
            double_bonds = 0
            for i in range(len(chain)-1):
                bond = mol.GetBondBetweenAtoms(chain[i], chain[i+1])
                if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                    double_bonds += 1
            # Our criteria: at least 12 carbons and at least 2 C=C bonds.
            if chain_length >= 12 and double_bonds >= 2:
                # If multiple candidate chains exist, choose the one with the greatest chain length,
                # or if equal, with more double bonds.
                if (chain_length > best_chain_info[0]) or (chain_length == best_chain_info[0] and double_bonds > best_chain_info[1]):
                    best_chain = chain
                    best_chain_info = (chain_length, double_bonds)
    
    # If no candidate chain met the criteria, report failure.
    if best_chain is None:
        # Optionally, provide some details why the candidates failed.
        reasons = []
        for anchor in candidate_indices:
            anchor_atom = mol.GetAtomWithIdx(anchor)
            for start in [nbr.GetIdx() for nbr in anchor_atom.GetNeighbors() if nbr.GetAtomicNum()==6 and nbr.GetIdx() in allowed]:
                chain = [anchor, start]
                current = start
                parent = anchor
                while True:
                    ext = [n for n in allowed.get(current, []) if n != parent]
                    if len(ext) == 1:
                        next_atom = ext[0]
                        chain.append(next_atom)
                        parent, current = current, next_atom
                    else:
                        break
                chain_length = len(chain)
                double_bonds = 0
                for i in range(len(chain)-1):
                    bond = mol.GetBondBetweenAtoms(chain[i], chain[i+1])
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        double_bonds += 1
                reasons.append(f"Chain length {chain_length} with {double_bonds} double bonds")
        return False, "No fatty acyl chain meets criteria: " + " | ".join(reasons)
    
    return True, (f"Contains fatty acyl chain with {best_chain_info[0]} carbons and "
                  f"{best_chain_info[1]} C=C bonds, meeting criteria for an essential fatty acid.")

# (Optional) Testing on sample SMILES strings.
if __name__ == '__main__':
    test_smiles = [
        "CC\\C=C/C\\C=C/C\\C=C/CCCCCC(O)=O",  # all-cis-7,10,13-hexadecatrienoic acid
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCC\\C=C/C\\C=C/C\\C=C/CCCCC",  # a phospholipid with a fatty acyl chain
        "OC(=O)CC\\C=C\\C\\C=C\\C\\C=C\\CCCCCCCC"  # example that might fail criteria
    ]
    for smi in test_smiles:
        result, reason = is_essential_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")