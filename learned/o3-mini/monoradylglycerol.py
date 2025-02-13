"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: Monoradylglycerol
Definition:
   Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position.
   
Our new approach:
 1. Parse the SMILES, add explicit hydrogens and sanitize.
 2. Find candidate glycerol backbones by scanning for a central sp3 carbon with two sp3 neighbors that are not connected to each other.
 3. For the three candidate glycerol carbons, require that each has exactly one oxygen substituent (i.e. oxygen neighbor not part of the candidate backbone).
 4. For each oxygen, determine if it is free (-OH, with at least one hydrogen) or substituted (no hydrogen).
 5. Exactly two free oxygens and one substituted oxygen must be found.
 6. For the substituted oxygen, follow the non-backbone branch and ensure there is a sufficiently long carbon chain (>=6 carbons) to indicate a fatty substituent.
 7. If a candidate passes these criteria, return True and a message.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines whether a molecule is a monoradylglycerol.
    A monoradylglycerol is defined as a glycerol (three contiguous sp3 carbons each bearing one oxygen substituent)
    in which exactly one of the oxygens is substituted by a fatty (acyl, alkyl or alk-1-enyl) chain, while the other 
    two remain as free hydroxyl (-OH) groups.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a monoradylglycerol, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so free hydroxyl groups show hydrogens explicitly.
    mol = Chem.AddHs(mol)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"
    
    # Helper: recursively compute the longest carbon chain length starting from an atom.
    def longest_chain_length(atom, visited):
        # Only count carbon atoms.
        if atom.GetAtomicNum() != 6:
            return 0
        max_length = 1
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                length = 1 + longest_chain_length(nbr, visited.copy())
                if length > max_length:
                    max_length = length
        return max_length

    # Gather all sp3 carbons.
    sp3_carbons = [atom for atom in mol.GetAtoms() 
                   if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3]
    
    # We will try to find a candidate glycerol backbone defined as:
    # A central sp3 carbon with two sp3 neighbors (the two terminal carbons) such that the two terminals are not connected.
    candidates_checked = []
    for center in sp3_carbons:
        # Get sp3 carbon neighbors of center.
        terminal_neighbors = [nbr for nbr in center.GetNeighbors() 
                              if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3]
        if len(terminal_neighbors) < 2:
            continue
        # Consider all pairs of distinct terminal neighbors.
        from itertools import combinations
        for term1, term2 in combinations(terminal_neighbors, 2):
            # In glycerol the two terminal carbons should not be directly bonded.
            if term1.GetIdx() in [a.GetIdx() for a in term2.GetNeighbors()]:
                continue
            backbone_idxs = sorted([center.GetIdx(), term1.GetIdx(), term2.GetIdx()])
            # Avoid rechecking the same backbone.
            if backbone_idxs in candidates_checked:
                continue
            candidates_checked.append(backbone_idxs)
            
            # Now, for each of the three backbone carbons, get oxygen substituents (neighbors not in the backbone).
            backbone_atoms = [center, term1, term2]
            oxygen_substituents = {}
            valid_backbone = True
            for atom in backbone_atoms:
                oxygens = [nbr for nbr in atom.GetNeighbors() 
                           if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [a.GetIdx() for a in backbone_atoms]]
                # For a clean glycerol, there should be exactly one oxygen substituent on each carbon.
                if len(oxygens) != 1:
                    valid_backbone = False
                    break
                oxygen_substituents[atom.GetIdx()] = oxygens[0]
            if not valid_backbone:
                continue
            
            # Now classify each oxygen as free (-OH) if it carries at least one hydrogen or substituted (acyl/ether) if not.
            free_count = 0
            substituted_data = None  # (parent carbon atom, oxygen) for the substituted one.
            for parent in backbone_atoms:
                ox = oxygen_substituents[parent.GetIdx()]
                # Using GetTotalNumHs() based on explicit hydrogens.
                if ox.GetTotalNumHs() > 0:
                    free_count += 1
                else:
                    # Assume this is the substituted oxygen.
                    substituted_data = (parent, ox)
            # Exactly two free hydroxyls and one substituted oxygen must be found.
            if free_count != 2 or substituted_data is None:
                continue
            
            # Check that the substituted oxygen is attached to a fatty chain.
            parent, sub_ox = substituted_data
            # Look for neighbors of sub_ox that are not the glycerol backbone.
            fatty_candidates = [nbr for nbr in sub_ox.GetNeighbors() 
                                if nbr.GetIdx() not in [a.GetIdx() for a in backbone_atoms]]
            if not fatty_candidates:
                continue
            fatty_found = False
            for candidate in fatty_candidates:
                if candidate.GetAtomicNum() == 6:  # a carbon must be part of the fatty substituent.
                    chain_length = longest_chain_length(candidate, set())
                    if chain_length >= 6:
                        fatty_found = True
                        break
            if not fatty_found:
                continue
            
            # If all conditions are met then we have found a monoradylglycerol.
            return True, "Found glycerol backbone with two free hydroxyls and one fatty substituent (monoradylglycerol)."
    
    # If no candidate backbone is found that qualifies:
    return False, "Could not identify a glycerol backbone with exactly one substituted oxygen representing a fatty chain."

# Example usage and testing:
if __name__ == "__main__":
    test_examples = [
        # Examples expected to be monoradylglycerols.
        "O(C(=O)CCCCCCCCC/C=C/CCCCCC)CC(O)CO",  # MG(18:1(11E)/0:0/0:0)[rac]
        "O(C(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)C(CO)CO",  # MG(0:0/22:5(...)/0:0)
        "C(\\C[C@H]1[C@H](CC([C@@H]1/C=C/[C@H](CCCCC)O)=O)O)=C\\CCCC(=O)OC(CO)CO",  # prostaglandin D2 2-glyceryl ester
        "CCCCCCCCCCCCCC\\C=C\\OC[C@@H](O)CO",  # 1-O-[(E)-hexadecen-1-yl]-sn-glycerol
        "CCCCCCCCC(=O)OC(CO)CO",  # 2-decanoylglycerol 
        "O(C(=O)CCCCCCCCCCCCCCCCCCCCC)C(CO)CO",  # MG(0:0/22:0/0:0)
        "CCCC(=O)OCC(O)CO",  # monobutyrin, should fail because fatty branch is too short
    ]
    for s in test_examples:
        result, reason = is_monoradylglycerol(s)
        print("SMILES:", s)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 50)