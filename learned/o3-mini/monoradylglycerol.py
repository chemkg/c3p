"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: Monoradylglycerol
Definition:
   Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent
   at an unspecified position.
Our revised approach:
 1. Parse the SMILES and add explicit hydrogens.
 2. Look over all sp3 carbon atoms to find a potential glycerol backbone.
    The candidate backbone is defined by a central sp3 carbon with at least two sp3 carbon neighbors (the terminals)
    which themselves are not directly connected.
 3. For each of these three carbons (backbone atoms), collect any oxygen neighbors (outside the backbone).
 4. If a backbone carbon has no oxygen neighbor, this candidate fails.
 5. When there is more than one oxygen neighbor, choose one—if one oxygen lacks hydrogens (i.e. is substituted) then choose that one.
 6. In the candidate backbone the classification criteria is that exactly two of the chosen oxygens are free hydroxyl groups
    (i.e. have at least one hydrogen) and exactly one is substituted (no hydrogen).
 7. For the substituted oxygen, follow its non-backbone neighbor (which should be carbon) and compute
    the longest chain length. If that chain has at least 6 carbons, this is taken as a fatty substituent.
 8. If any candidate passes, return True with a message. Otherwise, return False.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from itertools import combinations

def is_monoradylglycerol(smiles: str):
    """
    Determines whether a molecule is a monoradylglycerol.
    A monoradylglycerol is defined as a glycerol backbone in which exactly one oxygen has been substituted
    by a fatty acyl/alkyl/alk-1-enyl chain, while the other two oxygens remain as free hydroxyl (-OH) groups.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a monoradylglycerol, False otherwise.
        str: Explanation message for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that free -OH groups have H shown.
    mol = Chem.AddHs(mol)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"
    
    # Helper: compute the longest carbon chain length from a starting carbon atom,
    # counting only carbon atoms. Uses DFS with a visited set.
    def longest_chain_length(atom, visited):
        if atom.GetAtomicNum() != 6:
            return 0
        max_length = 1
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                chain_len = 1 + longest_chain_length(nbr, visited.copy())
                if chain_len > max_length:
                    max_length = chain_len
        return max_length

    # Gather all sp3 carbons (atomic number 6 and sp3 hybridized).
    sp3_carbons = [atom for atom in mol.GetAtoms() 
                   if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3]
    
    # Try to find a candidate glycerol backbone.
    # In glycerol the backbone consists of three contiguous carbons,
    # with a central carbon having at least two sp3 carbon neighbors that are not connected to each other.
    candidates_checked = []
    for center in sp3_carbons:
        # Find sp3 carbon neighbors (potential terminal carbons)
        terminals = [nbr for nbr in center.GetNeighbors() 
                     if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3]
        if len(terminals) < 2:
            continue
        for term1, term2 in combinations(terminals, 2):
            # In a glycerol backbone the two terminal carbons should not be directly bonded.
            if term1.GetIdx() in [a.GetIdx() for a in term2.GetNeighbors()]:
                continue
            # Assemble candidate backbone atoms (unique set of three indices).
            backbone_atoms = [center, term1, term2]
            backbone_ids = sorted(atom.GetIdx() for atom in backbone_atoms)
            if backbone_ids in candidates_checked:
                continue
            candidates_checked.append(backbone_ids)
            
            # For each backbone carbon, we expect at least one oxygen neighbor that is not in the backbone.
            chosen_oxygens = {}  # mapping backbone atom idx -> chosen oxygen atom
            candidate_valid = True
            for b_atom in backbone_atoms:
                o_neighbors = [nbr for nbr in b_atom.GetNeighbors() 
                               if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in backbone_ids]
                if not o_neighbors:
                    candidate_valid = False
                    break
                # If there is more than one oxygen attached, try to choose one that appears substituted (no hydrogen).
                chosen = None
                for o in o_neighbors:
                    # GetTotalNumHs() will include explicit H from the AddHs step.
                    if o.GetTotalNumHs() == 0:
                        chosen = o
                        break
                if chosen is None:
                    # Otherwise choose the first oxygen (which is free -OH)
                    chosen = o_neighbors[0]
                chosen_oxygens[b_atom.GetIdx()] = chosen
            if not candidate_valid:
                continue
            
            # Now count free hydroxyl groups and substituted oxygens in the backbone.
            free_count = 0
            substituted_info = None  # tuple (parent_carbon, oxygen)
            for b_atom in backbone_atoms:
                o_atom = chosen_oxygens[b_atom.GetIdx()]
                if o_atom.GetTotalNumHs() > 0:
                    free_count += 1
                else:
                    # If more than one is substituted then candidate fails.
                    if substituted_info is not None:
                        substituted_info = None
                        break
                    substituted_info = (b_atom, o_atom)
            # We need exactly 2 free hydroxyls and 1 substituted oxygen.
            if free_count != 2 or substituted_info is None:
                continue
            
            # Check that the substituted oxygen is attached to a valid fatty chain.
            parent_atom, sub_ox = substituted_info
            # Look for neighbor(s) of sub_ox that are not the backbone (should lead to the fatty chain)
            fatty_neighbors = [nbr for nbr in sub_ox.GetNeighbors() if nbr.GetIdx() not in backbone_ids]
            if not fatty_neighbors:
                continue
            fatty_found = False
            for fatty in fatty_neighbors:
                if fatty.GetAtomicNum() == 6:  # looking for a carbon branch
                    chain_len = longest_chain_length(fatty, set())
                    if chain_len >= 6:
                        fatty_found = True
                        break
            if not fatty_found:
                continue
            
            # If all conditions are met then we have found a monoradylglycerol.
            return True, "Found glycerol backbone with two free hydroxyls and one fatty substituent (monoradylglycerol)."
    
    return False, "Could not identify a glycerol backbone with exactly one substituted oxygen representing a fatty chain."

# Example usage and testing:
if __name__ == "__main__":
    test_examples = [
        "O(C(=O)CCCCCCCCC/C=C/CCCCCC)CC(O)CO",  # MG(18:1(11E)/0:0/0:0)[rac]
        "O(C(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)C(CO)CO",  # MG(0:0/22:5(...)/0:0)
        "C(\\C[C@H]1[C@H](CC([C@@H]1/C=C/[C@H](CCCCC)O)=O)O)=C\\CCCC(=O)OC(CO)CO",  # prostaglandin D2 2-glyceryl ester
        "CCCCCCCCCCCCCC\\C=C\\OC[C@@H](O)CO",  # 1-O-[(E)-hexadecen-1-yl]-sn-glycerol
        "CCCCCCCCC(=O)OC(CO)CO",  # 2-decanoylglycerol
        "O(C(=O)CCCCCCCCCCCCCCCCCCCCC)C(CO)CO",  # MG(0:0/22:0/0:0)
        "CCCC(=O)OCC(O)CO",  # monobutyrin (should fail due to short fatty chain)
        "OC(CO)COC(=O)CCCCCCC/C=C/C/C=C/C/C=C/CC",  # additional structure test
    ]
    for s in test_examples:
        result, reason = is_monoradylglycerol(s)
        print("SMILES:", s)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 50)