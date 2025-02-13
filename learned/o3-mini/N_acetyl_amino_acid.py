"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: N-acetyl-amino acid
Definition: An N-acyl-amino acid that has acetyl as the acyl group.
The molecule must have an N-acetyl group – that is, a nitrogen bearing an acetyl fragment 
(where the acyl carbon has a double-bonded oxygen and a connected methyl group) – and an amino 
acid backbone where that nitrogen is directly attached to an alpha-carbon that in turn is bonded 
to a carboxyl carbon (a carbon bonded to at least one double-bonded oxygen and one single-bonded oxygen).
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl amino acid based on its SMILES string.

    The method checks:
    1. That at least one nitrogen atom in the molecule carries an acetyl group (i.e. is attached to 
       a carbonyl carbon that is further attached to a methyl group).
    2. That the same nitrogen is also bonded to an alpha-carbon which is further connected to a carboxyl carbon.
       The carboxyl carbon is identified by having at least one oxygen bonded via a double bond (carbonyl) 
       and at least one oxygen bonded via a single bond.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is an N-acetyl amino acid, False otherwise.
        str: Reason for the classification result.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # Skip non-nitrogen atoms
        
        n_atom = atom  # candidate acetylated nitrogen
        acyl_found = None

        # Check neighbors of the nitrogen for an acetyl (acyl) carbon.
        # The acyl carbon should be bonded to n_atom,
        # have at least one double-bonded oxygen (carbonyl),
        # and be attached to a methyl group.
        for nbr in n_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue  # Look for a carbon neighbor
            acyl_carbon = nbr
            # Skip if this carbon is not bonded to the nitrogen (it is, by construction)
            double_oxygen = False
            methyl_found = False
            # Loop over neighbors of the acyl carbon (except the nitrogen)
            for nbr2 in acyl_carbon.GetNeighbors():
                if nbr2.GetIdx() == n_atom.GetIdx():
                    continue
                # Check for a double-bonded oxygen (the carbonyl oxygen)
                bond = mol.GetBondBetweenAtoms(acyl_carbon.GetIdx(), nbr2.GetIdx())
                if nbr2.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    double_oxygen = True
                # Check for a methyl group: a carbon with only one heavy (non-hydrogen) neighbor.
                if nbr2.GetAtomicNum() == 6:
                    # We use GetDegree() which counts the number of heavy-atom neighbors.
                    # A methyl carbon is usually only attached to one heavy atom (the acyl carbon)
                    if nbr2.GetDegree() == 1:
                        methyl_found = True
            if double_oxygen and methyl_found:
                acyl_found = acyl_carbon
                break  # Found an acetyl group on n_atom

        if acyl_found is None:
            # This nitrogen does not have an acetyl group; check the next nitrogen.
            continue

        # Now verify the amino acid backbone connectivity.
        # The nitrogen should be attached to an alpha-carbon (other than the acyl carbon).
        alpha_carbon = None
        for nbr in n_atom.GetNeighbors():
            if nbr.GetIdx() == acyl_found.GetIdx():
                continue  # skip the acetyl group carbon
            if nbr.GetAtomicNum() == 6:
                alpha_carbon = nbr
                break
        if alpha_carbon is None:
            continue  # No alpha-carbon attached to this nitrogen

        # Check that the alpha carbon has an attached carboxyl carbon.
        # In an amino acid, the alpha carbon is bound to a carboxyl carbon.
        carboxyl_found = False
        for nbr in alpha_carbon.GetNeighbors():
            # Skip the nitrogen that we already have.
            if nbr.GetIdx() == n_atom.GetIdx():
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            carboxyl_carbon = nbr
            # Evaluate carboxyl carbon: It should have at least one oxygen attached via a double bond 
            # and at least one oxygen attached via a single bond.
            double_count = 0
            single_count = 0
            for oxy in carboxyl_carbon.GetNeighbors():
                if oxy.GetAtomicNum() != 8:
                    continue
                bond = mol.GetBondBetweenAtoms(carboxyl_carbon.GetIdx(), oxy.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    double_count += 1
                elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    single_count += 1
            if double_count >= 1 and single_count >= 1:
                carboxyl_found = True
                break
        if carboxyl_found:
            return True, ("Found valid N-acetyl amino acid: N-atom has an acetyl group (with carbonyl and methyl) "
                          "and is bonded to an alpha-carbon which in turn is connected to a carboxyl group.")
    # If no nitrogen satisfied the connectivity requirements, return False.
    return False, "No N-acetylated nitrogen with an amino acid backbone (alpha-carbon with attached carboxyl group) was found"


# Example usage
if __name__ == "__main__":
    # We test with a few SMILES examples of N-acetyl amino acids.
    test_examples = {
        "N-acetyl-L-aspartic acid": "CC(=O)N[C@@H](CC(O)=O)C(O)=O",
        "N-acetylglycine": "CC(=O)NCC(O)=O",
        "N-acetylhistidine": "C(=O)(C(CC=1NC=NC1)NC(=O)C)O"
    }
    for name, smiles in test_examples.items():
        result, reason = is_N_acetyl_amino_acid(smiles)
        print(f"Name: {name}\nSMILES: {smiles}\nResult: {result}\nReason: {reason}\n{'-'*60}")