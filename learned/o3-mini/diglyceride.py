"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: diglyceride 
Definition: A glyceride that is glycerol in which any two of the hydroxy groups have been acylated.
That is, two positions carry an ester-linked oxygen (i.e. O–C(=O)–) and one position is "free"
(i.e. a hydroxyl or substituted by an alkyl group, but not esterified).
"""
from rdkit import Chem

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    The algorithm first searches for a glycerol-like backbone (three sp3 carbons connected in a chain)
    where leaving groups (substituents not in the backbone) occur exactly once per carbon.
    Then, for each of these substituents, it checks whether the oxygen is esterified (i.e. attached to a carbonyl).
    If exactly 2 substituents are ester oxygen (acylated) and the remaining one is free, the molecule is a diglyceride.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if it is a diglyceride, False otherwise.
        str: Reason for the classification.
    """
    # Parse molecule and add explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Helper function: regardless of stereochemistry, check if an oxygen substituent is part of an ester:
    # It is acylated if its neighbor (other than the candidate glycerol carbon) is a carbon that carries
    # at least one double bond to an oxygen.
    def is_acyl_oxygen(o_atom, glycerol_carbon_idx):
        for nbr in o_atom.GetNeighbors():
            if nbr.GetIdx() == glycerol_carbon_idx:
                continue
            if nbr.GetAtomicNum() == 6:  # carbon
                for bond in nbr.GetBonds():
                    # Look for a double bond to oxygen.
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    # We now search for candidate glycerol backbones.
    # A glycerol backbone should have three sp3 carbons connected in a chain
    # (i.e. one of them—the middle—is connected to the other two, which are not directly connected).
    num_atoms = mol.GetNumAtoms()
    checked_backbones = set()  # to avoid duplicate candidates
    
    for atom_i in mol.GetAtoms():
        if atom_i.GetAtomicNum() != 6 or atom_i.GetHybridization().name != "SP3":
            continue
        i = atom_i.GetIdx()
        # Look at neighbor carbons as a potential middle of the chain.
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
                # The chain should be linear. If the terminal carbons (i and k) are directly bonded, skip.
                if mol.GetBondBetweenAtoms(i, k) is not None:
                    continue
                backbone = tuple(sorted([i, j, k]))
                if backbone in checked_backbones:
                    continue
                checked_backbones.add(backbone)
                
                # Identify the middle atom: it should be the one connected to the other two.
                backbone_set = set(backbone)
                middle_atom = None
                terminals = []
                for idx in backbone:
                    atom = mol.GetAtomWithIdx(idx)
                    n_in_backbone = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() in backbone_set]
                    if len(n_in_backbone) == 2:
                        middle_atom = atom
                    elif len(n_in_backbone) == 1:
                        terminals.append(atom)
                if middle_atom is None or len(terminals) != 2:
                    continue  # not a well-defined chain
                
                # Now check each candidate glycerol carbon.
                # In a “pure” glycerol backbone, each carbon should have exactly one substituent (i.e.
                # an atom not in the backbone). If there are extra connections (for instance, if another heavy atom
                # like phosphorus is attached) we skip this candidate.
                acylated_positions = 0
                free_positions = 0
                valid_backbone = True
                reasons = []
                for c_idx in backbone:
                    atom = mol.GetAtomWithIdx(c_idx)
                    # Get non-backbone heavy neighbors (atomic number > 1 and not in backbone)
                    subs = [nbr for nbr in atom.GetNeighbors() 
                            if nbr.GetIdx() not in backbone and nbr.GetAtomicNum() > 1]
                    if len(subs) != 1:
                        valid_backbone = False
                        reasons.append(f"Atom {c_idx} does not have exactly one heavy substituent (found {len(subs)})")
                        break
                    sub = subs[0]
                    # If the substituent is phosphorus or connected to phosphorus, then skip candidate.
                    if sub.GetAtomicNum() == 15:
                        valid_backbone = False
                        reasons.append(f"Substituent on atom {c_idx} is phosphorus")
                        break
                    # If the substituent is oxygen, decide if it is acylated.
                    if sub.GetAtomicNum() == 8:
                        if is_acyl_oxygen(sub, c_idx):
                            acylated_positions += 1
                        else:
                            free_positions += 1
                    else:
                        # If not oxygen, we treat it as free (for example an alkyl group).
                        free_positions += 1
                # We require exactly 2 acylated and 1 free positions.
                if valid_backbone:
                    if acylated_positions == 2 and free_positions == 1:
                        reason_found = ("Found glycerol backbone (three contiguous sp3 carbons) with "
                                        "two acylated substituents (ester-linked oxygens) and one free group.")
                        return True, reason_found
    # If no candidate passes our filters, return False with message.
    return False, "No valid glycerol backbone with exactly two acylated positions and one free position found; not a diglyceride."

# Example usage (testing with one of the provided DG structures)
if __name__ == '__main__':
    test_smiles = "O(C(=O)CCCCCCCCCCCCC(C)C)[C@H](COC(=O)CCCCCCCCCCCCC(C)C)CO"  # DG(i-16:0/i-16:0/0:0)
    result, reason = is_diglyceride(test_smiles)
    print(result, reason)