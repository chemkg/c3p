"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: A derivative of the dimethylisoalloxazine skeleton (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione)
with a substituent on the 10 position (i.e. a flavin).
Note: This SMARTS-based approach is an approximation.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin (a derivative of the dimethylisoalloxazine skeleton
    with a substituent on the 10 position) based on its SMILES string.

    The method works by:
      1. Parsing the SMILES.
      2. Checking for the presence of an isoalloxazine core that bears the 7,8-dimethyl pattern.
         The SMARTS used is:
           [n:10]1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O
         Here the atom tagged [n:10] is taken as the N at position 10.
      3. For any match the degree (number of bonds) of the atom corresponding to N10 is compared
         to the degree in the bare core (query). In the query the N10 has only the two bonds 
         needed for the ring closure. In a true flavin there is an extra substituent attached to N10.
    
    Args:
        smiles (str): SMILES string for the molecule.

    Returns:
        bool: True if the molecule matches the flavin criteria, False otherwise.
        str: A reason describing the classification decision.
    """
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define SMARTS for the dimethylisoalloxazine core.
    # This query requires that the core contains two methyl groups on the benzene ring 
    # (i.e. the 7,8-dimethyl pattern). We label the nitrogen at the 10 position using [n:10].
    core_smarts = "[n:10]1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in the SMARTS pattern for the flavin core."
    
    # Check if the molecule contains the isoalloxazine core.
    # Note: GetSubstructMatches returns a tuple of tuples, each with the atom indices in mol that match the query.
    matches = mol.GetSubstructMatches(core_query)
    if not matches:
        return False, "Isoalloxazine core with 7,8-dimethyl pattern not found."
    
    # Identify which atom in the query has the mapping number of 10.
    n10_query_idx = None
    for atom in core_query.GetAtoms():
        if atom.HasProp("molAtomMapNumber") and atom.GetProp("molAtomMapNumber") == "10":
            n10_query_idx = atom.GetIdx()
            break
    if n10_query_idx is None:
        return False, "The query pattern does not contain an atom mapped with number 10."
    
    # Get the degree (number of bonds) in the bare isoalloxazine core for the N10 atom.
    core_n10_degree = core_query.GetAtomWithIdx(n10_query_idx).GetDegree()
    # In our query the tagged N10 is only connected to 2 atoms (the ring bonds).
    
    # Now, check the actual substitution at the corresponding N10 in the molecule.
    # We will loop over all matches (if the core occurs more than once, one may qualify as flavin).
    for match in matches:
        # match is a tuple giving the indices in mol corresponding to the query atoms.
        n10_mol_idx = match[n10_query_idx]
        n10_atom = mol.GetAtomWithIdx(n10_mol_idx)
        # Compare the number of neighbors. If the actual degree is greater than in the query,
        # it indicates that N10 carries an extra substituent.
        if n10_atom.GetDegree() > core_n10_degree:
            return True, "Molecule contains the expected dimethylisoalloxazine core with a substituent on N10."
    
    return False, "Isoalloxazine core found but no substituent attached at the 10 position."

# Example usage (uncomment for testing):
# smiles_list = [
#     "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C",  # lumiflavin
#     "[H]C(=O)n1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O",  # 7,8-dimethylisoalloxazine-10-carbaldehyde
#     "CC1=CC2=C(C=C1C)N(C3=NC(=O)NC(=O)C3=N2)C[C@H](O)[C@H](O)[C@H](O)CO"  # riboflavin derivative
# ]
# for sm in smiles_list:
#     result, reason = is_flavin(sm)
#     print(sm, "->", result, "-", reason)