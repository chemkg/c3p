"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: ether lipid
Definition: A lipid similar in structure to a glycerolipid but in which one or more of the carbon atoms on glycerol 
            is bonded to an alkyl chain via an ether linkage instead of the usual ester linkage.
            
This program refines previous attempts by broadening the glycerol backbone detection (allowing for acylation)
and further requiring that an ether oxygen connects to a long (aliphatic) chain.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    It first checks if a glycerol-like backbone is present (using two alternative SMARTS patterns)
    and then verifies that at least one alkyl chain is attached to the glycerol via a non-ester (ether) oxygen.
    Furthermore, one of the carbon neighbors of the ether oxygen must extend into a long aliphatic chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an ether lipid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization failed: " + str(e)
    
    # --- 1. Detect a glycerol-like backbone.
    #
    # Instead of a single pattern, we allow either a free triol backbone or one where one hydroxyl is replaced by a carbonyl:
    glycerol_patterns = [
        Chem.MolFromSmarts("OCC(O)CO"),    # free glycerol backbone (all hydroxyls)
        Chem.MolFromSmarts("OCC(=O)CO")     # one substitution (e.g. ester at sn-2)
    ]
    glycerol_found = False
    for patt in glycerol_patterns:
        if mol.HasSubstructMatch(patt):
            glycerol_found = True
            break
    if not glycerol_found:
        return False, "No glycerol backbone detected (using patterns OCC(O)CO or OCC(=O)CO)"
    
    # --- 2. Define a helper to decide whether an oxygen is in an ester group.
    def is_ester_oxygen(o_atom):
        # An ester oxygen is typically attached to a carbon that is double-bonded to an oxygen.
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon
                for bond in nbr.GetBonds():
                    # bond.GetBondTypeAsDouble() returns 2.0 for a double bond.
                    if bond.GetBondTypeAsDouble() == 2.0:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    # --- 3. Helper function to estimate the maximum aliphatic chain length starting from a given carbon.
    # We perform a simple DFS along carbon-carbon bonds (ignoring rings).
    def max_chain_length_from(atom, visited):
        length = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                # To favor aliphatic chains, we skip atoms in rings.
                if nbr.IsInRing():
                    continue
                visited.add(nbr.GetIdx())
                branch_length = 1 + max_chain_length_from(nbr, visited)
                if branch_length > length:
                    length = branch_length
                visited.remove(nbr.GetIdx())
        return length

    # --- 4. Search for a candidate ether oxygen.
    # We look for an oxygen atom that is:
    #   - connected to exactly two heavy atoms,
    #   - NOT carrying any explicit hydrogen,
    #   - not part of an ester linkage,
    #   - and at least one of its carbon neighbors gives rise to a long aliphatic chain.
    ether_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # oxygen
            # Check that it is not an -OH (should have no explicit hydrogen; note: implicit hydrogens not counted in neighbors)
            if atom.GetTotalNumHs() > 0:
                continue
            
            heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
            if len(heavy_neighbors) != 2:
                continue

            # Skip if the oxygen appears to be part of an ester linkage.
            if is_ester_oxygen(atom):
                continue

            # Both neighbors must be carbons
            if not all(nbr.GetAtomicNum() == 6 for nbr in heavy_neighbors):
                continue

            # For at least one neighbor, the attached aliphatic chain should be "long" (here 6 or more carbons).
            chain_ok = False
            for c in heavy_neighbors:
                # We allow the neighbor itself to be counted as the start of the chain.
                chain_length = 1 + max_chain_length_from(c, {c.GetIdx()})
                if chain_length >= 6:
                    chain_ok = True
                    break

            if chain_ok:
                ether_found = True
                break

    if not ether_found:
        return False, "No non-ester (ether) linkage found attached to glycerol that leads to a long alkyl chain"
    
    return True, "Molecule contains a glycerol backbone with at least one alkyl chain attached via an ether linkage"

# Example usage:
if __name__ == "__main__":
    # Use one of the provided examples:
    test_smiles = "CCCCCCCCCCCCCC\\C=C\\OC[C@@H](O)CO"  # 1-[(E)-hexadecen-1-yl]-sn-glycerol
    result, reason = is_ether_lipid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)