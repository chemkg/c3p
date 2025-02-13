"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine

Definition:
An alkyl,acyl-sn-glycero-3-phosphocholine in which unspecified alkyl and acyl groups are located 
at positions 1 and 2 respectively. Such molecules comprise three parts:
    • An O-alkyl (ether) chain (position 1)
    • An O-acyl (ester) chain (position 2)
    • A phosphocholine headgroup (attached via phosphate at position 3)

This program uses RDKit substructure searches as a heuristic to identify these moieties.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    
    The function checks that:
      - The molecule contains a phosphocholine headgroup, using the SMARTS pattern:
            P(=O)([O-])OCC[N+](C)(C)C
      - The molecule contains at least one ester (acyl) moiety (i.e. an oxygen attached to C=O)
      - The molecule contains an O-alkyl (ether) linkage that is not part of a phosphate or ester.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule passes the tests and is deemed a member of this lipid class.
        str: A reason string for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for phosphocholine headgroup.
    # This pattern matches a phosphate group with a phosphocholine fragment.
    phos_pattern = Chem.MolFromSmarts("P(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phos_pattern):
        return False, "Missing phosphocholine headgroup (pattern P(=O)([O-])OCC[N+](C)(C)C not found)"

    # 2. Look for an ester group corresponding to an acyl chain.
    # The ester group is indicated by an oxygen attached to a carbonyl: OC(=O)
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "Missing ester (acyl) group (pattern OC(=O) not found)"

    # 3. Look for an O-alkyl (ether) bond that is not part of a phosphate or an ester.
    # We do this by iterating over oxygen atoms in the molecule.
    alkoxy_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:  # only consider oxygen atoms
            continue
        # Skip oxygens that are attached to phosphorus (they are part of the phosphate headgroup)
        if any(neigh.GetAtomicNum() == 15 for neigh in atom.GetNeighbors()):
            continue

        # Check the neighbors of this oxygen.
        # We want an oxygen that is connected to at least one carbon that is not involved in a carbonyl.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon neighbor
                # Check whether this carbon is a carbonyl carbon.
                is_carbonyl = False
                for bond in nbr.GetBonds():
                    # If the bond is a double bond and the other atom is oxygen, treat it as carbonyl.
                    if bond.GetBondTypeAsDouble() == 2:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            is_carbonyl = True
                            break
                if not is_carbonyl:
                    alkoxy_found = True
                    break
        if alkoxy_found:
            break
    if not alkoxy_found:
        return False, "Missing O-alkyl (ether) group for the alkyl chain"

    # Optional additional checks:
    # You might add further tests (e.g., verifying chain lengths, number of heteroatoms, etc.)
    # For now, if all three chemical features are present, classify as the lipid type.
    
    return True, ("Contains phosphocholine headgroup, acyl ester linkage, and O-alkyl ether substituent "
                   "consistent with 2-acyl-1-alkyl-sn-glycero-3-phosphocholine.")

# Example usage:
if __name__ == "__main__":
    # one example from the provided list:
    test_smiles = "CCCCCCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCC\\C=C/CCCCCC"
    result, reason = is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(test_smiles)
    print(f"Test result: {result}\nReason: {reason}")