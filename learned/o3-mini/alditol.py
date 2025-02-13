"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: alditol 
Definition: A carbohydrate that is an acyclic polyol having the general formula 
HOCH2[CH(OH)]nCH2OH (formally derivable from an aldose by reduction of the carbonyl group).
This program attempts to find an acyclic carbon chain whose pattern meets:
  - Terminal carbons are CH2OH (exactly one heavy neighbor from the chain plus one oxygen)
  - Interior carbons are CH(OH) (two chain neighbors plus one oxygen)
and the chain does not belong to any ring.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule can be classified as an alditol based on its SMILES string.
    The approach:
      1) Rejects molecules containing a carbonyl group.
      2) Searches for a contiguous acyclic chain of carbons decorated with exactly one –OH per carbon,
         with the terminal carbons being CH2OH and the interior ones CH(OH).
         (The chain must not be part of any ring.)
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule matches an alditol pattern, False otherwise.
        str: A text explanation of the decision.
    """
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for consistent neighbor counting.
    mol = Chem.AddHs(mol)
    
    # Rule out molecules that have a carbonyl group (C=O) because that indicates an unreduced aldose.
    carbonyl_query = Chem.MolFromSmarts("[#6]=[OX1]")
    if mol.HasSubstructMatch(carbonyl_query):
        return False, "Contains a carbonyl group, so not a fully reduced (alditol) structure"
    
    # Try finding a chain (with different possible lengths)
    # We'll consider chains from 3 up to 12 carbons.
    for n in range(3, 13):
        # Build a SMARTS string for the chain:
        # Terminal: CH2OH is represented as "CO"
        # Interior: CH(OH) is represented as "C(O)"
        # The expression below builds a linear chain: "CO" + "C(O)"*(n-2) + "CO"
        pattern_smarts = "CO" + "C(O)" * (n - 2) + "CO"
        query = Chem.MolFromSmarts(pattern_smarts)
        if query is None:
            continue  # in case pattern is invalid
        
        # Search for substructure matches without using chirality (to be more flexible).
        matches = mol.GetSubstructMatches(query, useChirality=False)
        for match in matches:
            # 'match' is in the same order as the pattern, so it represents a candidate chain in order.
            # We now verify that:
            #   a) None of the chain carbons is in a ring.
            #   b) Each chain carbon has exactly the expected heavy neighbors:
            #      - For a terminal carbon (first and last): 1 chain neighbor and 1 non-chain oxygen (total heavy count ==2).
            #      - For an interior carbon: 2 chain neighbors and 1 non-chain oxygen (total heavy count ==3).
            valid_chain = True
            chain_set = set(match)
            # Loop over the chain atoms in the order provided by the match.
            for i, atom_idx in enumerate(match):
                atom = mol.GetAtomWithIdx(atom_idx)
                # (a) Reject if the atom is in a ring.
                if atom.IsInRing():
                    valid_chain = False
                    break
                # Get heavy (non-H) neighbors.
                heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                # Determine how many neighbors in the chain are expected.
                if i == 0 or i == (len(match) - 1):
                    expected_chain_count = 1  # terminal: only one adjacent chain carbon
                else:
                    expected_chain_count = 2  # interior: two adjacent chain carbons
                # Identify which of the heavy neighbors are part of the chain.
                chain_neighbors = [nbr.GetIdx() for nbr in heavy_neighbors if nbr.GetIdx() in chain_set]
                # The extra neighbor should be the hydroxyl oxygen.
                non_chain_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetIdx() not in chain_set]
                
                # Check: the number of chain neighbors should equal the expected count.
                if len(chain_neighbors) != expected_chain_count:
                    valid_chain = False
                    break
                # And there should be exactly one non-chain heavy neighbor and it must be oxygen.
                if len(non_chain_neighbors) != 1 or non_chain_neighbors[0].GetAtomicNum() != 8:
                    valid_chain = False
                    break
            if valid_chain:
                return True, f"Contains an acyclic polyol chain matching HOCH2[CH(OH)]^{n-2}CH2OH pattern"
    
    return False, "No acyclic polyol chain matching the required HOCH2[CH(OH)]nCH2OH pattern was found"

# For simple testing when run as a script:
if __name__ == "__main__":
    # A few SMILES examples: erythritol and glycerol.
    test_smiles = [
        "OC[C@H](O)[C@H](O)CO",   # erythritol (should be True)
        "OCC(O)CO",               # glycerol (3 carbons but typically not considered an alditol from an aldose reduction)
        "OC[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)CO"  # D-iditol (should be True)
    ]
    for smi in test_smiles:
        result, reason = is_alditol(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")