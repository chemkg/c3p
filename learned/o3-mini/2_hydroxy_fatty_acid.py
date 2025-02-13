"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: Any fatty acid that has a hydroxy functional group in the alpha- or 2-position,
i.e. a 2-hydroxy fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is defined as a fatty acid containing a carboxylic acid group
    and an -OH substituent on the carbon alpha (adjacent) to that carboxyl group.
    
    The procedure is:
      1. Parse the SMILES and add explicit hydrogens.
      2. Check that the molecule has at least one carboxylic acid group.
      3. For each carboxyl group found, retrieve the carbon atom and examine its (non-oxygen)
         neighbors – i.e. candidates for the alpha carbon. For each such candidate, check
         if it has any oxygen atom attached that carries an explicit hydrogen (i.e. an -OH).
      4. Also enforce that the molecule is acyclic (common for simple fatty acids) and that
         it has a sufficient number of carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 2-hydroxy fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES to get an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that our search can detect -OH groups.
    molH = Chem.AddHs(mol)
    
    # Many fatty acids are linear (acyclic). Reject molecules with rings.
    if molH.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s), not a linear fatty acid"
    
    # Count carbon atoms (we require at minimum 4 carbons; this can be adjusted)
    carbon_count = sum(1 for atom in molH.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4:
        return False, "Too few carbons to be considered a fatty acid"
    
    # Define SMARTS for a carboxylic acid group.
    carboxyl_smarts = "C(=O)[O;H1,O-]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = molH.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) == 0:
        return False, "Missing a carboxylic acid group"
    
    # For at least one carboxyl group, check the alpha carbon for an -OH substituent.
    # Each match from carboxyl_pattern gives atom indices, where the first is the carbon.
    alpha_found = False
    for match in carboxyl_matches:
        carboxyl_c_idx = match[0]  # carbon atom in the carboxyl group
        carboxyl_c = molH.GetAtomWithIdx(carboxyl_c_idx)
        # Iterate over each neighbor of the carboxyl carbon.
        for neighbor in carboxyl_c.GetNeighbors():
            # We want a carbon that is not one of the oxygens in the carboxyl.
            if neighbor.GetAtomicNum() != 6:
                continue
            # For the candidate alpha carbon, check if it has an -OH group
            for sub_neigh in neighbor.GetNeighbors():
                # We are looking for an oxygen with at least one hydrogen attached.
                if sub_neigh.GetAtomicNum() == 8 and sub_neigh.GetTotalNumHs() > 0:
                    alpha_found = True
                    break
            if alpha_found:
                break
        if alpha_found:
            break

    if not alpha_found:
        return False, "Missing alpha (2-) hydroxy group adjacent to a carboxyl group"
    
    # Optional: Check that the overall molecular weight is in a fatty acid range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(molH)
    if mol_wt < 100:  # Most fatty acids are over about 100 Da.
        return False, "Molecular weight too low for a fatty acid"
    
    return True, "Molecule is a 2-hydroxy fatty acid with an alpha-hydroxy group adjacent to the carboxyl group"

# Example usage:
# test_smiles = "CCCCCCCC\\C=C/CCCCCCCCCCCCC(O)C(O)=O"  # 2-hydroxynervonic acid
# result, reason = is_2_hydroxy_fatty_acid(test_smiles)
# print(result, reason)