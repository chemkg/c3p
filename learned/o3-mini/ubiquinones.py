"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: Ubiquinones – any benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone.
Ubiquinones usually carry a polyprenyl side chain at position 6.
"""

from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule belongs to the class of ubiquinones based on its SMILES string.
    The classification is based on the presence of a quinone ring derived from 2,3-dimethoxy-5-methylbenzoquinone,
    along with (optionally) an isoprenoid/polyprenyl side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a ubiquinone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS query matching the ubiquinone core --
    # This pattern represents a 6-membered ring with two carbonyl groups 
    # (indicative of the quinone) and with two methoxy groups and one methyl group.
    # Note: The SMARTS "O=C1C(OC)=C(OC)C(=O)C(C)=C1" is our approximation to the 2,3-dimethoxy-5-methylbenzoquinone skeleton.
    ubiquinone_core_smarts = "O=C1C(OC)=C(OC)C(=O)C(C)=C1"
    core_query = Chem.MolFromSmarts(ubiquinone_core_smarts)
    if core_query is None:
        # In case the SMARTS pattern itself does not parse, return None.
        return None, None

    if not mol.HasSubstructMatch(core_query):
        return False, "Does not contain the 2,3-dimethoxy-5-methylbenzoquinone core"
    
    # Optionally, search for an isoprene-like unit indicative of a polyprenyl side chain.
    # This SMARTS looks for a fragment roughly matching an isoprenoid unit ("C/C=C(C)").
    prenyl_smarts = "C/C=C(C)"
    prenyl_query = Chem.MolFromSmarts(prenyl_smarts)
    n_prenyl = len(mol.GetSubstructMatches(prenyl_query)) if prenyl_query is not None else 0

    # Build a reason string. Even if no prenyl unit is detected,
    # the ubiquinone core is the defining feature, and short-chain homologues do exist.
    reason = "Contains 2,3-dimethoxy-5-methylbenzoquinone core"
    if n_prenyl > 0:
        reason += f" with {n_prenyl} isoprenoid unit(s) detected"
    else:
        reason += " but no isoprenoid side chain detected (this likely indicates a short-chain ubiquinone homolog)"
    
    return True, reason