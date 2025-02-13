"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: N-acetyl-amino acid
Definition: An N-acyl-amino acid that has acetyl as the acyl group.
This version first tries to identify an amino acid backbone by locating a free carboxyl group 
(with a carbonyl and an –OH that is not negatively charged) and its adjacent alpha carbon. 
Then every nitrogen in the molecule is checked to see if it carries an acetyl group – that is, 
if it is directly bonded to a carbon that is part of a CH3-C(=O)- group. Finally, for each acetylated 
nitrogen the shortest-bond distance to any alpha carbon is computed, and if that distance is within 
a (somewhat relaxed) threshold (here up to 6 bonds) we classify the molecule as an N-acetyl-amino acid.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl amino acid based on its SMILES string.
    
    The algorithm follows these steps:
    (1) Parse the SMILES string.
    (2) Identify a free carboxyl group. A carboxyl carbon is defined as a carbon 
        that has one oxygen attached by a double bond and one oxygen attached by 
        a single bond where that oxygen has degree 1 and no formal charge.
    (3) Identify candidate alpha carbons – carbons that are directly bonded 
        to the carboxyl carbon and that also have at least one neighboring nitrogen.
    (4) For every nitrogen atom in the molecule, test if it is acetylated. 
        We consider a nitrogen acetylated if it is bonded to a carbon for which:
            - At least one neighboring oxygen is double-bonded (the carbonyl O)
            - At least one other neighbor is a methyl group (a carbon with only one heavy-atom neighbor)
    (5) For every acetylated nitrogen found, compute the shortest path (in bonds)
        between that nitrogen and any candidate alpha carbon. If the distance is less than
        or equal to 6 bonds, we consider the acetyl group to be “backbone-connected” and 
        classify the molecule as an N‑acetyl amino acid.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule qualifies as an N‑acetyl amino acid, False otherwise.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Helper to test if a carbon atom is a "free" carboxyl carbon.
    def is_carboxyl_carbon(carbon_atom):
        # Must be carbon.
        if carbon_atom.GetAtomicNum() != 6:
            return False
        dblO_found = False
        oh_found = False
        # Loop over neighbors to check for carbonyl oxygen and –OH oxygen.
        for neigh in carbon_atom.GetNeighbors():
            if neigh.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), neigh.GetIdx())
            if bond is None:
                continue
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                dblO_found = True
            elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                # In a free carboxyl group, the –OH oxygen should have degree==1 
                # and be not charged.
                if neigh.GetDegree() == 1 and neigh.GetFormalCharge() == 0:
                    oh_found = True
        return dblO_found and oh_found

    # Step (2): find candidate carboxyl carbons and corresponding alpha carbons.
    alpha_carbon_idxs = set()
    for atom in mol.GetAtoms():
        # Look for carboxyl carbon
        if atom.GetAtomicNum() != 6:
            continue
        if not is_carboxyl_carbon(atom):
            continue
        carboxyl_idx = atom.GetIdx()
        # For each neighbor of this carboxyl carbon, if it is carbon,
        # and if that neighbor is attached to any nitrogen, consider it as an alpha carbon.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                # Check if any neighbor (other than the carboxyl C) is a nitrogen.
                for sub_nbr in nbr.GetNeighbors():
                    if sub_nbr.GetIdx() == carboxyl_idx:
                        continue
                    if sub_nbr.GetAtomicNum() == 7:
                        alpha_carbon_idxs.add(nbr.GetIdx())
                        break

    if not alpha_carbon_idxs:
        return False, "No valid amino acid backbone (free carboxyl group with adjacent alpha carbon) found"

    # Helper to test if a given nitrogen is acetylated.
    # The idea: the nitrogen should have a carbon neighbor (the acyl carbon) that in turn
    # has (i) at least one oxygen attached by a double bond and (ii) a methyl group neighbor.
    def is_acetylated_nitrogen(n_atom):
        for nbr in n_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue  # must be carbon
            acyl_c = nbr
            carbonyl_found = False
            methyl_found = False
            for subnbr in acyl_c.GetNeighbors():
                # Skip going back to the nitrogen
                if subnbr.GetIdx() == n_atom.GetIdx():
                    continue
                if subnbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(acyl_c.GetIdx(), subnbr.GetIdx())
                    if bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        carbonyl_found = True
                elif subnbr.GetAtomicNum() == 6:
                    # Check if this is likely a methyl group.
                    # We assume that the methyl carbon is only bonded to the acyl carbon.
                    if subnbr.GetDegree() == 1:
                        methyl_found = True
            if carbonyl_found and methyl_found:
                return True
        return False

    # Step (4): Collect indices of all atoms that are acetylated nitrogen.
    acetylated_n_idxs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        if is_acetylated_nitrogen(atom):
            acetylated_n_idxs.append(atom.GetIdx())
    if not acetylated_n_idxs:
        return False, "No nitrogen bearing an acetyl group found"

    # Step (5): For each acetylated nitrogen, check connectivity to the amino acid backbone.
    # We require that the shortest path (in bonds) from the acetylated N to at least one
    # candidate alpha carbon is at most 6 bonds.
    for n_idx in acetylated_n_idxs:
        for a_idx in alpha_carbon_idxs:
            path = rdmolops.GetShortestPath(mol, n_idx, a_idx)
            # path is a tuple of atom indices; length-1 gives number of bonds
            if 1 <= (len(path) - 1) <= 6:
                return True, ("Found valid N-acetyl amino acid: An acetylated nitrogen (bearing a CH3-C(=O)- fragment) "
                              "is connected to an amino acid backbone (free carboxyl group attached to an alpha carbon) "
                              f"via a path of {len(path)-1} bonds.")
    return False, "No acetylated nitrogen is sufficiently connected to an amino acid backbone (free carboxyl group adjacent to an alpha carbon) found"

# Example usage:
if __name__ == "__main__":
    test_smiles = {
        "grixazone A": "[H]C(=O)c1ccc2oc3cc(=O)c(N)c(SC[C@H](NC(C)=O)C(O)=O)c3nc2c1",
        "N-acetylhistidine": "C(=O)(C(CC=1NC=NC1)NC(=O)C)O",
        "N-acetyl-L-aspartic acid": "CC(=O)N[C@@H](CC(O)=O)C(O)=O",
        "bisorcic": "C(CCNC(=O)C)[C@@H](C(O)=O)NC(=O)C",
        "N-acetyl-D-methionine": "CSCC[C@@H](NC(C)=O)C(O)=O",
        "N-acetyl-beta-alanine": "CC(=O)NCCC(O)=O",
        "N(6)-acetyl-L-lysine": "CC(=O)NCCCC[C@H](N)C(O)=O",
        "N(5)-acetyl-L-ornithine": "CC(=O)NCCC[C@H](N)C(O)=O",
        "N-acetyl-L-proline": "CC(=O)N1CCC[C@H]1C(O)=O",
        "N-acetyl-L-phenylalaninate (false positive test)": "CC(=O)N[C@@H](Cc1ccccc1)C([O-])=O",
    }
    for name, sm in test_smiles.items():
        result, reason = is_N_acetyl_amino_acid(sm)
        print(f"Name: {name}\nSMILES: {sm}\nResult: {result}\nReason: {reason}\n{'-'*60}")