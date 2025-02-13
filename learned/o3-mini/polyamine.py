"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: Polyamine
Definition: Any organic amino compound that contains two or more free amino groups.
A free amino group is defined here as a nitrogen atom (atomic number 7) that has at least 
one hydrogen attached (after making H atoms explicit) and is not directly bound to a carbonyl
group (i.e. is not part of an amide bond). This helps to avoid falsely classifying peptides.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is any organic compound that contains at least two free amino groups.
    (A free amino group here is a nitrogen with at least one hydrogen attached and not forming an
    amide bond with an adjacent carbonyl carbon.)
    
    Additionally, if the molecule appears to be a peptide (i.e. contains one or more amide bonds)
    and is sufficiently large, it will not be classified as a polyamine.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polyamine, False otherwise.
        str: Reason for the classification result.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check organic: require at least one carbon atom
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule does not contain any carbon atoms, thus not organic"
    
    # Work with explicit hydrogens so that hydrogen counts are accessible
    mol = Chem.AddHs(mol)
    
    # Heuristic: if the molecule contains amide bonds then it might be a peptide or polyamide.
    # We use the simple SMARTS pattern for an amide bond: C(=O)N
    amide_query = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_query)
    # Also get molecular weight as computed by RDKit
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if amide_matches and mol_wt > 300:
        return False, "Contains amide bonds and high molecular weight, likely a peptide or polyamide"
    
    # Helper function: decide whether a nitrogen atom is in an amide bond.
    def is_amide_nitrogen(n_atom, mol_obj):
        # We assume n_atom is a nitrogen.
        # Check if any neighbor carbon is bound (via a single bond) to this nitrogen and 
        # in turn has a double bond to an oxygen. That carbonyl is indicative of an amide bond.
        for nbr in n_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon neighbor
                bond = mol_obj.GetBondBetweenAtoms(n_atom.GetIdx(), nbr.GetIdx())
                if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Check if this carbon is bonded to any oxygen by a double bond.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetAtomicNum() == 8:
                            bond2 = mol_obj.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                            if bond2 and bond2.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                return True
        return False
    
    # Count free amino groups: nitrogen atoms with at least one (explicit) hydrogen that are not
    # part of an amide bond.
    free_amino_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        # Get total number of attached hydrogens (explicit only)
        h_count = atom.GetTotalNumHs()
        if h_count < 1:
            continue
        # Exclude those nitrogens that are part of an amide bond
        if is_amide_nitrogen(atom, mol):
            continue
        free_amino_count += 1
    
    if free_amino_count < 2:
        return False, f"Contains {free_amino_count} free amino group(s), need at least 2 for polyamine"
    else:
        return True, f"Contains {free_amino_count} free amino groups, satisfying polyamine criteria"

# Example usage (can be removed in production):
if __name__ == "__main__":
    # Test with trimethylenediamine as a simple polyamine example
    test_smiles = "NCCCN"  # trimethylenediamine
    result, reason = is_polyamine(test_smiles)
    print(result, reason)