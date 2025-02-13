"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: N-acetyl-amino acid
Definition: An N-acyl-amino acid that has acetyl as the acyl group.
This program checks whether a molecule (given by its SMILES string) has
1) at least one nitrogen that is acetylated – that is, the nitrogen is bonded
   to a carbon that carries a carbonyl oxygen and a methyl group.
2) such a nitrogen is “in the amino acid backbone”, meaning that there is a short (2–3 bond) path
   from that N to a carboxyl group. Here the carboxyl group is defined as a carbon
   that is double-bonded to an oxygen and that is singly bonded to an oxygen that is not further substituted by a carbon.
The use of a distance check (via a shortest-path search) rather than simply insisting on a direct bond 
(from nitrogen to an alpha–carbon) allows for both classical (α‑) amino acids as well as β‑amino acids.
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl amino acid based on its SMILES string.
    
    Steps:
      1. Parse the SMILES.
      2. For every nitrogen atom in the molecule, check if it bears an acetyl group.
         An acetyl group here is defined as a carbon neighbor that has at least one
         double-bonded oxygen (i.e. a carbonyl) and also a methyl group (a carbon neighbor
         with only one heavy-atom connection).
      3. For each N-acetylated nitrogen, search for a carboxyl group in the molecule that
         is “backbone-connected” – that is, the shortest path (in bonds) from the N to that carboxyl carbon is 2 or 3 bonds.
         The carboxyl carbon is defined to have:
             - at least one oxygen attached by a double bond
             - at least one oxygen attached by a single bond that is not substituted by carbons 
               (i.e. the oxygen’s only heavy neighbor is the carboxyl carbon).
    If such a nitrogen is found, the function returns (True, <explanation>).
    Otherwise, it returns (False, <reason>).
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if molecule qualifies, False otherwise.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # helper function to test if a carbon is part of a free carboxyl group.
    def is_carboxyl_carbon(carbon_atom):
        # Ensure the atom is carbon.
        if carbon_atom.GetAtomicNum() != 6:
            return False
        has_carbonyl = False
        has_OH = False
        # Loop over neighbors looking for oxygens.
        for neigh in carbon_atom.GetNeighbors():
            if neigh.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), neigh.GetIdx())
            if bond is None:
                continue
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                has_carbonyl = True
            elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                # In a free acid, the single-bonded oxygen should be –OH.
                # We assume that if the oxygen is attached only to the carboxyl carbon (degree==1)
                # then it is an –OH rather than an –OR.
                if neigh.GetDegree() == 1:
                    has_OH = True
        return has_carbonyl and has_OH

    # Iterate over nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # skip non-nitrogen atoms
        n_atom = atom
        acetyl_found = False
        acyl_neighbor = None  # the carbon of the acetyl group
        # Look through neighbors of the nitrogen for an acetyl pattern:
        for nbr in n_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue  # must be carbon
            acyl_candidate = nbr
            double_oxygen = False
            methyl_found = False
            # Check neighbors of the acyl candidate (skip back the nitrogen)
            for subnbr in acyl_candidate.GetNeighbors():
                if subnbr.GetIdx() == n_atom.GetIdx():
                    continue
                if subnbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(acyl_candidate.GetIdx(), subnbr.GetIdx())
                    if bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        double_oxygen = True
                elif subnbr.GetAtomicNum() == 6:
                    # A methyl carbon typically is only connected to this acyl candidate.
                    if subnbr.GetDegree() == 1:
                        methyl_found = True
            if double_oxygen and methyl_found:
                acetyl_found = True
                acyl_neighbor = acyl_candidate
                break
        if not acetyl_found:
            continue  # this nitrogen is not acetylated; check next nitrogen

        # Now check for an amino acid backbone.
        # Instead of demanding the N is directly attached to the alpha-carbon,
        # we allow for either standard (α-) amino acids (path length=2) or β‑amino acids (path length=3).
        # We look for any carboxyl carbon (c, as defined) for which the shortest bond path from
        # the N (our acetylated N) is 2 or 3.
        carboxyl_found = False
        for atom2 in mol.GetAtoms():
            # Look at carbon atoms that can be carboxyl groups.
            if atom2.GetAtomicNum() != 6:
                continue
            if not is_carboxyl_carbon(atom2):
                continue
            # Find shortest path (list of atom indices) from our N to this carboxyl carbon.
            path = rdmolops.GetShortestPath(mol, n_atom.GetIdx(), atom2.GetIdx())
            # The number of bonds is len(path)-1; for backbone connectivity we allow 2 or 3 bonds.
            if 2 <= (len(path) - 1) <= 3:
                carboxyl_found = True
                break
        if carboxyl_found:
            return True, ("Found valid N-acetyl amino acid: The nitrogen bears an acetyl group "
                          "and a carboxyl group is connected via a short path (2–3 bonds) from that nitrogen.")
    return False, "No N-acetylated nitrogen with a suitable amino acid backbone (free carboxyl group connected within 2–3 bonds) was found"


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
        "N-acetyl-L-proline": "CC(=O)N1CCC[C@H]1C(O)=O",
        "N-acetyl-L-phenylalaninate": "CC(=O)N[C@@H](Cc1ccccc1)C([O-])=O",  # expected false positive
    }
    for name, sm in test_smiles.items():
        result, reason = is_N_acetyl_amino_acid(sm)
        print(f"Name: {name}\nSMILES: {sm}\nResult: {result}\nReason: {reason}\n{'-'*60}")