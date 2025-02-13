"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: Galactosylceramide – any cerebroside in which the monosaccharide head group is galactose.
In this version, rather than using a SMARTS with the quantifier {10,} (which failed), we 
examine each amide carbonyl to find the acyl chain and then walk the chain to count consecutive CH2 groups.
A chain of at least 10 CH2 groups (on one side of an amide) is considered typical for a ceramide.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is defined as a glycosphingolipid in which the ceramide (a sphingoid base
    linked via an amide bond to a long-chain fatty acid) bears a single galactose.

    This function applies the following filters:
      1. The molecule must contain a galactose ring (alpha or beta form) as the sugar head group.
      2. The molecule must contain an amide bond.
      3. Among any amide bonds, at least one must have an acyl chain that includes at least
         10 consecutive CH2 groups.
      4. Crude overall filters on molecular weight and carbon count are applied.

    Args:
       smiles (str): SMILES string of the molecule.

    Returns:
       bool: True if the molecule is classified as a galactosylceramide, False otherwise.
       str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Helper function: determine if an atom is a CH2 group.
    def is_ch2(atom):
        # For our purposes, a CH2 is a carbon (not aromatic) with exactly two attached hydrogens.
        if atom.GetAtomicNum() != 6:
            return False
        # Use GetTotalNumHs which takes into account implicit hydrogens.
        # Also require that the atom is not aromatic.
        return (atom.GetTotalNumHs() == 2) and (not atom.GetIsAromatic())

    # Helper function: recursively walk the chain of CH2 groups.
    def longest_chain(atom, prev):
        max_length = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == prev.GetIdx():
                continue
            if is_ch2(nbr):
                # Continue the chain from this neighbor.
                length = 1 + longest_chain(nbr, atom)
                if length > max_length:
                    max_length = length
        return max_length

    # 1. Check for galactose head group.
    # Define two substructures for galactose (alpha and beta forms).
    beta_gal = Chem.MolFromSmiles("CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O")
    alpha_gal = Chem.MolFromSmiles("CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O")
    has_gal = False
    if beta_gal and mol.HasSubstructMatch(beta_gal):
        has_gal = True
    elif alpha_gal and mol.HasSubstructMatch(alpha_gal):
        has_gal = True
    if not has_gal:
        return False, "No galactose sugar head group found"

    # 2. Check for an amide bond that is part of a ceramide.
    # Instead of relying solely on a SMARTS with a quantifier, we now iterate over atoms to find an amide carbonyl.
    fatty_acyl_found = False
    # Loop over atoms to look for a carbonyl carbon in an amide group.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Look for a double bond to oxygen.
        oxygen_bond = None
        nitrogen_found = False
        for bond in atom.GetBonds():
            # Check for a double bond to an oxygen.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    oxygen_bond = nbr
            # Also check for a single bond to a nitrogen (the amide nitrogen).
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 7:
                    nitrogen_found = True
        # If we have both a carbonyl (oxygen_bond) and an amide nitrogen, consider this atom.
        if oxygen_bond is not None and nitrogen_found:
            # Identify the acyl chain—one of the neighbors that is neither the oxygen (from C=O)
            # nor the nitrogen (from the amide bond).
            acyl_neighbor = None
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() in [7, 8]:
                    continue
                acyl_neighbor = nbr
                break
            # If an acyl chain candidate exists and is a CH2 group, follow the chain.
            if acyl_neighbor and is_ch2(acyl_neighbor):
                # Count this CH2 and recursively count the following consecutive CH2 units.
                chain_length = 1 + longest_chain(acyl_neighbor, atom)
                if chain_length >= 10:
                    fatty_acyl_found = True
                    break

    if not fatty_acyl_found:
        return False, "No long fatty acyl chain (>=10 consecutive CH2 groups) detected on the amide – not a typical ceramide"

    # 3. Apply crude overall filters.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low to be a typical glycosphingolipid or ceramide derivative"
    # Count carbons.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbon atoms to represent a typical ceramide aliphatic chain"

    return True, ("Contains a galactose sugar head group and an amide-linked ceramide backbone "
                  "with a long fatty acyl chain (galactosylceramide)")
                    
# Example tests (uncomment to run):
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_galactosylceramide(test_smiles)
# print(result, reason)