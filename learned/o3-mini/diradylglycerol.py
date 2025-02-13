"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: Diradylglycerol
Definition: Any lipid that is glycerol bearing two substituent groups – either acyl, alkyl, or alk-1-enyl – 
at any two of the three possible positions.
Heuristic:
  1. Parse the molecule from the SMILES.
  2. Look for a contiguous 3-carbon chain (candidate glycerol backbone).
  3. For each of the three carbons, get the oxygen(s) attached that are not part of the backbone.
     In natural glycerol (C3H8O3) there is one –OH per carbon. In a diradylglycerol two of these –OH’s are substituted 
     with acyl/alkyl/alk-1-enyl groups (and thus the oxygen will usually not carry an explicit hydrogen), 
     while the remaining position should have a free hydroxyl ([OX2H]).
  4. If we find a candidate glycerol backbone with exactly 3 oxygen neighbors – two of which are substituted 
     (i.e. no hydrogens) and one free (has at least one hydrogen) – we return True.
  
If any of these steps fail the heuristic, we return False along with a reason.
"""

from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string according to the heuristic:
      - It must include a glycerol backbone (three connected carbons)
      - Exactly two of the positions are substituted (with acyl, alkyl or alk-1-enyl groups) whilst
        the third retains a free hydroxyl group.
        
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a diradylglycerol; False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # We will search for a glycerol backbone: a chain of 3 carbon atoms (not necessarily labeled)
    # We then inspect the oxygens attached to these carbons (ignoring carbons that are in the backbone).
    # For a true glycerol unit there should be 3 oxygen substituents (one at each carbon).
    # In a diradylglycerol exactly two of these oxygens are substituted (no hydrogen) while the remaining one is free (-OH).
    
    # Get all carbon atoms indices.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "C"]
    
    # Flag to indicate if a glycerol with proper substitution was found.
    found_glycerol = False
    reason = ""
    
    # Iterate over all atoms and try to identify a 3-carbon chain.
    # For each carbon, look for two neighbors that are carbons to possible form a chain.
    n_atoms = mol.GetNumAtoms()
    for a in range(n_atoms):
        atom_a = mol.GetAtomWithIdx(a)
        if atom_a.GetSymbol() != "C":
            continue
        # Look for a neighbor b that is carbon.
        for b in [nbr.GetIdx() for nbr in atom_a.GetNeighbors() if nbr.GetSymbol() == "C"]:
            atom_b = mol.GetAtomWithIdx(b)
            # b should be connected to a and then to c (a chain of three carbons)
            for c in [nbr.GetIdx() for nbr in atom_b.GetNeighbors() if nbr.GetSymbol() == "C" and nbr.GetIdx() not in (a,)]:
                # Ensure the chain is in order: a-b-c. (order does not matter)
                chain_indices = {a, b, c}
                # For each of these three carbons, get oxygen neighbors that are not in the chain.
                oxy_details = []  # list of tuples: (atom_index, is_free_OH: bool)
                valid_chain = True
                for idx in (a, b, c):
                    carbon_atom = mol.GetAtomWithIdx(idx)
                    # Get oxygen neighbors that are not part of the glycerol chain.
                    oxy_neighbors = [nbr for nbr in carbon_atom.GetNeighbors() if nbr.GetSymbol() == "O" and nbr.GetIdx() not in chain_indices]
                    # In natural glycerol, each carbon should have exactly one oxygen attached.
                    if len(oxy_neighbors) != 1:
                        valid_chain = False
                        break
                    # For this oxygen, determine if it is a free hydroxyl.
                    oxy_atom = oxy_neighbors[0]
                    # Using GetTotalNumHs() to decide if the oxygen carries any hydrogen(s).
                    numHs = oxy_atom.GetTotalNumHs()
                    # If at least one H, we treat it as a free hydroxyl; otherwise, it is substituted.
                    is_free = (numHs > 0)
                    oxy_details.append((oxy_atom.GetIdx(), is_free))
                if not valid_chain:
                    continue
                
                # Count substituted vs free hydroxyls in the candidate glycerol fragment.
                free_count = sum(1 for (_, is_free) in oxy_details if is_free)
                sub_count = sum(1 for (_, is_free) in oxy_details if not is_free)
                # For a diradylglycerol, we expect exactly 1 free hydroxyl and 2 substituted positions.
                if free_count == 1 and sub_count == 2:
                    found_glycerol = True
                    reason = ("Found glycerol backbone with two substituted positions (likely acyl/alkyl/alk-1-enyl) "
                              "and one free hydroxyl.")
                    return True, reason
                else:
                    reason = (f"Found glycerol backbone but substitution pattern is freeOH: {free_count}, "
                              f"substituted: {sub_count} (requires 1 free hydroxyl and 2 substituted).")
    
    if not found_glycerol:
        return False, ("No glycerol backbone with exactly two substituent groups (and one free hydroxyl) was identified. "
                       "The molecule may not be a diradylglycerol.")
    
    return False, "Unexpected end of classification routine."

# Example usage:
if __name__ == "__main__":
    # Test the function with one of the provided SMILES strings (DG(20:4(5Z,8Z,11Z,14Z)/18:2(9Z,12Z)/0:0))
    test_smiles = "C([C@@](COC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)[H])O"
    classification, reason = is_diradylglycerol(test_smiles)
    print("Classification:", classification)
    print("Reason:", reason)