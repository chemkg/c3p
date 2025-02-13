"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: diglyceride 
Definition: A glyceride that is glycerol in which any two of the hydroxy groups have been acylated.
That is, two positions carry an ester-linked oxygen (O–C(=O)–) and one position is free (–OH or substituted by an alkyl group).
"""

from rdkit import Chem

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    The algorithm first rejects molecules that contain phosphorus since those tend to be complex phospholipids.
    It then searches for candidate glycerol backbones (three contiguous sp3 carbons forming a chain with
    one heavy substituent per carbon). For each such candidate, the substituents are evaluated: an oxygen substituent 
    is marked as acylated if it is attached (besides the glycerol carbon) to a carbon that carries at least one double 
    bond to oxygen (i.e. a carbonyl group). Any substituent that is not oxygen is treated as free.
    A diglyceride should have exactly 2 acylated positions and 1 free position.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if it is classified as a diglyceride, False otherwise.
        str: Reason for the classification.
    """
    # Parse molecule; if parsing fails, return
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add hydrogens to aid in substructure analysis
    mol = Chem.AddHs(mol)
    
    # Quickly rule out molecules with phosphorus (common in phospholipids)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            return False, "Contains phosphorus, likely not a diglyceride"
    
    # Helper: check if an oxygen substituent is part of an ester bond.
    # It is acylated if it is attached (besides the bonded glycerol carbon) to a carbon that has at least one C=O bond.
    def is_acyl_oxygen(o_atom, glycerol_carbon_idx):
        for nbr in o_atom.GetNeighbors():
            if nbr.GetIdx() == glycerol_carbon_idx:
                continue
            if nbr.GetAtomicNum() == 6:  # candidate acyl carbon
                # Look for a double bond from this carbon to an oxygen
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    # Search for candidate glycerol backbones.
    # A glycerol backbone is considered to be three sp3 hybridized carbons connected in a chain,
    # with one of them (the middle) bound to both terminal carbons
    num_atoms = mol.GetNumAtoms()
    checked_backbones = set()  # avoid duplicate candidates
    
    for atom_i in mol.GetAtoms():
        if atom_i.GetAtomicNum() != 6 or atom_i.GetHybridization().name != "SP3":
            continue
        i = atom_i.GetIdx()
        # Consider neighbor carbons as potential middle carbon candidates.
        for atom_j in atom_i.GetNeighbors():
            if atom_j.GetAtomicNum() != 6 or atom_j.GetHybridization().name != "SP3":
                continue
            j = atom_j.GetIdx()
            for atom_k in atom_j.GetNeighbors():
                if atom_k.GetAtomicNum() != 6 or atom_k.GetHybridization().name != "SP3":
                    continue
                k = atom_k.GetIdx()
                if k == i:
                    continue
                # For a linear chain, ensure that terminals (i and k) are not directly bonded.
                if mol.GetBondBetweenAtoms(i, k) is not None:
                    continue
                backbone = tuple(sorted([i, j, k]))
                if backbone in checked_backbones:
                    continue
                checked_backbones.add(backbone)
                
                # Identify the middle and terminal carbons.
                backbone_set = set(backbone)
                middle_atom = None
                terminals = []
                for idx in backbone:
                    atom = mol.GetAtomWithIdx(idx)
                    nbrs_in_backbone = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() in backbone_set]
                    if len(nbrs_in_backbone) == 2:
                        middle_atom = atom
                    elif len(nbrs_in_backbone) == 1:
                        terminals.append(atom)
                if middle_atom is None or len(terminals) != 2:
                    continue  # not a valid three-carbon chain
                
                # Now, in a glycerol fragment (with no extra linkage), each backbone carbon should have exactly one heavy substituent outside of the backbone.
                acylated_positions = 0
                free_positions = 0
                candidate_valid = True
                reasons = []
                
                for c_idx in backbone:
                    atom = mol.GetAtomWithIdx(c_idx)
                    # Count substituents that are not in the backbone (ignoring hydrogens).
                    subs = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in backbone and nbr.GetAtomicNum() > 1]
                    if len(subs) != 1:
                        candidate_valid = False
                        reasons.append(f"Atom index {c_idx} does not have exactly one heavy substituent (found {len(subs)})")
                        break
                    sub = subs[0]
                    # If the substituent is phosphorus, skip candidate.
                    if sub.GetAtomicNum() == 15:
                        candidate_valid = False
                        reasons.append(f"Substituent on atom index {c_idx} is phosphorus")
                        break
                    # If the substituent is oxygen then determine if it is esterified (acylated) or free.
                    if sub.GetAtomicNum() == 8:
                        if is_acyl_oxygen(sub, c_idx):
                            acylated_positions += 1
                        else:
                            free_positions += 1
                    else:
                        # Non-oxygen substituents are treated as free (e.g. alkyl groups).
                        free_positions += 1
                # Candidate backbone must present exactly two acylated and one free position
                if candidate_valid:
                    if acylated_positions == 2 and free_positions == 1:
                        return True, ("Found glycerol backbone (three contiguous sp3 carbons) with "
                                      "two acylated substituents (ester-linked oxygens) and one free group.")
    return False, "No valid glycerol backbone with exactly two acylated positions and one free position found; not a diglyceride."

# Example usage:
if __name__ == '__main__':
    # Test with one of the provided DG structures:
    test_smiles = "O(C(=O)CCCCCCCCCCCCC(C)C)[C@H](COC(=O)CCCCCCCCCCCCC(C)C)CO"  # DG(i-16:0/i-16:0/0:0)
    result, reason = is_diglyceride(test_smiles)
    print(result, reason)