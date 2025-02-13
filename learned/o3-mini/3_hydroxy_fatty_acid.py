"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: 3-hydroxy fatty acids
Definition: Any fatty acid with a hydroxy functional group in the beta- (or 3-) position.
The molecule must contain a terminal carboxylic acid group (i.e. -C(=O)O) and a hydroxyl
(-OH) attached at the beta-carbon (i.e. two bonds away from the carboxyl carbon).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid must have a terminal carboxylic acid and a hydroxyl group
    on the beta (3-) position relative to the acid group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to detect -OH properly.
    mol = Chem.AddHs(mol)
    
    # Detect carboxylic acid group.
    # SMARTS: [CX3](=O)[OX2H1] matches a typical carboxylic acid carbon.
    acid_smarts = "[CX3](=O)[OX2H1]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # Try each carboxylic acid found - we expect a fatty acid to have one terminal acid.
    for match in acid_matches:
        # In the SMARTS, match[0] is the carboxyl carbon.
        acid_carbon_idx = match[0]
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
        
        # Find the neighbor carbon (alpha carbon) attached to the carboxyl carbon.
        # Exclude the oxygens (which are part of the acid group).
        alpha_candidates = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not alpha_candidates:
            continue  # Try next acid group if available.
        # Assuming linear fatty acid, take the first alpha candidate.
        alpha_atom = alpha_candidates[0]
        
        # Now, the beta (or 3-) carbon shall be a neighbor of the alpha carbon but not the acid carbon.
        beta_candidates = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetIdx() != acid_carbon_idx and nbr.GetAtomicNum() == 6]
        if not beta_candidates:
            continue
        
        # For each candidate beta carbon, check for an attached hydroxyl (-OH) group.
        for beta_atom in beta_candidates:
            for nbr in beta_atom.GetNeighbors():
                # We are looking for an oxygen attached via a single bond that carries at least one hydrogen.
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(beta_atom.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        # In the RDKit molecule with explicit hydrogens, the -OH oxygen should have H attached.
                        if nbr.GetTotalNumHs() >= 1:
                            return True, "Found a hydroxy (-OH) group attached at the beta (3-) position relative to a carboxylic acid group"
    
    return False, "No beta-hydroxy group found on the fatty acid chain"