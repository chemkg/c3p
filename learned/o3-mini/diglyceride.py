"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: diglyceride 
Definition: A glyceride that is glycerol in which any two of the hydroxy groups have been acylated.
That is, two of the positions have an ester group (an oxygen attached to a carbonyl carbon) while
the remaining position is free (OH or substituted with an alkyl group).
"""
from rdkit import Chem

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride is a glyceride in which exactly two of the three hydroxyl groups of glycerol
    are acylated (i.e. have been converted to esters). The remaining group is left free (OH or alkyl).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a diglyceride, False otherwise.
        str: Reason for the classification.
    """
    # Parse the molecule and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Helper function to check if an oxygen substituent is acylated.
    # We consider an oxygen (attached to a candidate glycerol carbon) to be acylated if 
    # it has a neighbor (other than the glycerol carbon) which is a carbon that is bonded via a double bond to an oxygen.
    def is_acyl_oxygen(o_atom, glycerol_carbon_idx):
        for nbr in o_atom.GetNeighbors():
            if nbr.GetIdx() == glycerol_carbon_idx:
                continue
            if nbr.GetAtomicNum() == 6:  # carbon
                # Check if this carbon has a double bond to oxygen.
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    # We now try to find a glycerol backbone candidate:
    # Look for a chain of 3 carbon atoms (not necessarily in a ring) where one carbon (the middle)
    # is connected to the two others.
    # We iterate over all triplets of carbon atoms that are connected in a linear manner.
    num_atoms = mol.GetNumAtoms()
    candidate_found = False
    reason_found = ""
    checked_backbones = set()  # to avoid duplicates
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        i = atom.GetIdx()
        # i is one carbon candidate; search its carbon neighbors for a potential middle atom.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            j = nbr.GetIdx()
            # Now for the neighbor j, search for a second carbon neighbor (other than i) that is also carbon.
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetAtomicNum() != 6:
                    continue
                k = nbr2.GetIdx()
                if k == i:
                    continue
                # Ensure that i and k are not directly bonded (since the three-carbons should form a chain, not a cycle)
                if mol.GetBondBetweenAtoms(i, k) is not None:
                    continue
                # Now, form a sorted tuple for uniqueness.
                backbone = tuple(sorted([i, j, k]))
                if backbone in checked_backbones:
                    continue
                checked_backbones.add(backbone)
                # Identify which atom is in the middle of the chain.
                # The middle should be the one that is connected to the other two.
                chain_atoms = {i, j, k}
                middle = None
                terminals = []
                for idx in backbone:
                    a = mol.GetAtomWithIdx(idx)
                    # Count how many neighbors of 'a' are in the chain.
                    count_in_chain = sum(1 for nbr in a.GetNeighbors() if nbr.GetIdx() in chain_atoms)
                    if count_in_chain == 2:
                        middle = a
                    elif count_in_chain == 1:
                        terminals.append(a)
                if middle is None or len(terminals) != 2:
                    # Not a proper chain of 3 carbons.
                    continue

                # For each of the three glycerol candidate carbons, check its substituents (neighbors not in the backbone).
                # In glycerol, we expect three substituents – one on each carbon. In diglyceride two should be acylated (ester oxygen)
                # and one should be free (likely a free hydroxyl or replaced by an alkyl group).
                acylated_positions = 0
                free_positions = 0
                # We will check each candidate carbon.
                for a in (mol.GetAtomWithIdx(i), mol.GetAtomWithIdx(j), mol.GetAtomWithIdx(k)):
                    # Find all neighbors not in the backbone that are oxygen atoms.
                    oxygen_neighbors = []
                    for nbr in a.GetNeighbors():
                        if nbr.GetIdx() not in backbone and nbr.GetAtomicNum() == 8:
                            oxygen_neighbors.append(nbr)
                    # If there is at least one oxygen neighbor, see if any qualifies as acylated.
                    found_acyl = False
                    for o in oxygen_neighbors:
                        if is_acyl_oxygen(o, a.GetIdx()):
                            found_acyl = True
                            break
                    if found_acyl:
                        acylated_positions += 1
                    else:
                        # If there is no oxygen substituent or the oxygen is not acylated,
                        # we treat this position as free.
                        free_positions += 1
                # We require exactly two acylated positions and one free position.
                if acylated_positions == 2 and free_positions == 1:
                    candidate_found = True
                    reason_found = ("Found glycerol backbone (three contiguous sp3 carbons) "
                                    "with two acylated substituents (ester-linked oxygens) and one free group.")
                    return True, reason_found

    if not candidate_found:
        return False, "No valid glycerol backbone with two acylated positions and one free position found; not a diglyceride."
    
    # If by some chance we break out, return False.
    return False, "Unknown failure in diglyceride classification."

# Example usage:
if __name__ == '__main__':
    # You can test with one of the provided SMILES strings for DG.
    test_smiles = "O(C(=O)CCCCCCCCCCCCC(C)C)[C@H](COC(=O)CCCCCCCCCCCCC(C)C)CO"
    result, reason = is_diglyceride(test_smiles)
    print(result, reason)