"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: A derivative of the dimethylisoalloxazine (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione)
skeleton with a substituent attached on the 10 position (i.e. a flavin).
This approach uses a two‐step strategy:
  1. A lenient SMARTS is used to find an isoalloxazine core.
  2. The candidate N atom (tagged in the query as N10 via [nH:10]) is inspected:
       • In the unsubstituted core it should have an attached hydrogen.
         In a flavin the extra substituent should replace that hydrogen.
       • We further require that the extra substituent attached is “carbon‐based”
         (as expected for a typical ribityl chain) to help discount direct adenine attachments.
Note: Being a SMARTS‐based heuristic, this is an approximation.
"""
from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin derivative (i.e. a derivative of the dimethylisoalloxazine
    skeleton with a substituent at the 10 position) based on its SMILES string.
    
    Our strategy is:
      1. Parse the SMILES.
      2. Look for the isoalloxazine core using a lenient SMARTS.
         (We do not force two methyl substituents on the benzene ring so as to catch potentially
         modified flavin derivatives.)
         The query is written as:
           "n1c2cc(cc2)nc2c1nc(=O)[nH:10]c2=O"
         where we tag the N that in the unsubstituted core would be N10.
      3. For any match, examine the candidate N (mapped as 10):
           - In the unsubstituted isoalloxazine the N would be drawn as [nH],
             so it should have one bound hydrogen.
           - In a proper flavin, an extra substituent is attached at this nitrogen so that
             the bound hydrogen is replaced (i.e. the N now has no H).
           - We also inspect the substituent(s) attached to N10 that are not part of the core.
             One of these should be a carbon (atomic number 6) representing the start
             of a ribityl or similar chain.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a flavin, False otherwise.
       str: A reason describing the decision.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a lenient SMARTS pattern to capture the isoalloxazine core.
    # We do not insist on two methyl groups on the benzene ring.
    # In an unsubstituted isoalloxazine the N at the 10-position is drawn as [nH].
    # We tag that [nH] with mapping number 10.
    core_smarts = "n1c2cc(cc2)nc2c1nc(=O)[nH:10]c2=O"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in the SMARTS pattern for the flavin core."
    
    # Find substructure matches of the isoalloxazine core in the molecule.
    matches = mol.GetSubstructMatches(core_query)
    if not matches:
        return False, "Isoalloxazine core not found."
    
    # For each match (if there are multiple, one flavin hit is enough)
    for match in matches:
        # First, identify the candidate atom corresponding to N10 in the query.
        # We find the query atom with mapping number "10" – note that our SMARTS only tags one atom.
        n10_query_idx = None
        for atom in core_query.GetAtoms():
            if atom.HasProp("molAtomMapNumber") and atom.GetProp("molAtomMapNumber") == "10":
                n10_query_idx = atom.GetIdx()
                break
        if n10_query_idx is None:
            return False, "The query pattern does not contain the expected atom mapping for N10."
    
        # Get the corresponding atom in the molecule.
        n10_mol_idx = match[n10_query_idx]
        n10_atom = mol.GetAtomWithIdx(n10_mol_idx)
    
        # In the bare (unsubstituted) core the N10 would be drawn as [nH] (i.e. with one hydrogen).
        # If it carries an extra substituent then it will have lost that hydrogen.
        # We check the total number of attached hydrogens (both explicit and implicit).
        n10_total_h = n10_atom.GetTotalNumHs()
        if n10_total_h > 0:
            # Still has an H, so no substituent is attached at N10.
            return False, "Isoalloxazine core found but no substituent attached at the 10 position."
    
        # Now, verify that at least one neighbor (outside the core) is a carbon.
        # (This is meant to check that an alkyl chain -- typically a ribityl chain -- is attached.)
        # First, obtain the set of atom indices that are part of the core match.
        core_atom_indices = set(match)
        valid_substituent = False
        for neighbor in n10_atom.GetNeighbors():
            if neighbor.GetIdx() not in core_atom_indices:
                # Check if the substituent is carbon-based.
                if neighbor.GetAtomicNum() == 6:
                    valid_substituent = True
                    break
        if not valid_substituent:
            return False, "Substituent at N10 is not carbon‐based as expected for a flavin."
    
        # We can also (if desired) add further criteria (for example, checking the overall size,
        # flexibility, or polar groups) to further discriminate out false positives.
    
        return True, "Molecule contains the expected isoalloxazine core with a substituent on N10."
    
    return False, "No valid flavin substructure match found."

# Example usage:
if __name__ == "__main__":
    # Test with a few examples:
    test_smiles = {
        "8-formyl-8-demethylriboflavin 5'-phosphate":
            "C(N1C=2C(=NC3=C1C=C(C(=C3)C)C(=O)[H])C(NC(N2)=O)=O)[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)O",
        "Riboflavol":
            "OC(CN1C2=C(N=C3C1=NC(=O)N(O)C3=O)C=C(C(=C2)C)C)C(O)C(O)CO",
        "riboflavin cyclic 4',5'-phosphate":
            "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C[C@H](O)[C@H](O)[C@H]3COP(O)(=O)O3)c2cc1C",
        "FAD (a correct flavin)":
            "Cc1cc2nc3c(nc(=O)[nH]c3=O)N(C[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)OP(O)(=O)OC[C@H]3O[C@H]([C@H](O)[C@@H]3O)n3cnc4c(N)ncnc34)c2cc1C",
        "Flavin adenine (false positive)":
            "P(OC[C@H]1O[C@@H]([N+]=2C3=NC=NC(N)=C3NC2)[C@H](O)[C@@H]1O)(OP(OCC(O)C(O)C(O)CN4C5=C(N=C6C4=NC(=O)NC6=O)C=C(C(=C5)C)C)(O)=O)(O)=O"
    }
    for name, smi in test_smiles.items():
        flag, reason = is_flavin(smi)
        print(f"{name}: {flag} – {reason}")