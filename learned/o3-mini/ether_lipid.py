"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: ether lipid
Definition: A lipid similar in structure to a glycerolipid but in which one or more of the carbon atoms on glycerol 
            is bonded to an alkyl chain via an ether (C–O–C) linkage instead of the usual ester (CO–O) linkage.
            
Our strategy:
1. Detect a glycerol “backbone.” In glycerol proper (HO–CH2–CHOH–CH2–OH) the backbone carbons can be
   pointed out by the SMARTS "OCC(O)CO". When one position is acylated we may see "OCC(=O)CO".
2. For each backbone match we extract the backbone carbon indices (for the free glycerol pattern, these
   are at positions 1, 2 and 3) and then check whether at least one of these carbons is attached to an oxygen
   that is not part of the backbone. That oxygen must not be part of an ester linkage (i.e. not adjacent to a C=O)
   and must lead into a long aliphatic chain (that we roughly assess via a DFS along C–C bonds excluding rings).

Note: This method is heuristic and may “miss” unusual representations or flag non-lipids that mimic the pattern.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    
    The algorithm:
      1. Parse the SMILES.
      2. Look for a glycerol-like backbone using two SMARTS patterns (free glycerol or monoacylated).
      3. For each glycerol match, take the backbone carbons (for a free glycerol match, the atoms
         at positions 1,2,3 are assumed to be the glycerol carbons). Then, inspect all neighbors of these carbons:
         if one of them is an oxygen (not already part of the backbone) and that oxygen is not “esterified”
         then check that from the oxygen the non-backbone side (a neighboring carbon) leads into a long
         aliphatic chain (using a DFS that ignores rings) with at least 6 carbons.
      4. If such a feature is found, the molecule is classified as an ether lipid.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an ether lipid, False otherwise.
        str: Explanation for the classification.
    """
    # --- Parse and sanitize the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization failed: " + str(e)
    
    # --- 1. Look for a glycerol-like backbone.
    # The free glycerol backbone: HO–CH2–CHOH–CH2–OH is represented as "OCC(O)CO" (5 atoms: O, C, C, C, O)
    # A monoacyl (esterified) backbone may appear as "OCC(=O)CO".
    glycerol_smarts_patterns = [
        ("OCC(O)CO", [1, 2, 3]),  # for free glycerol, backbone carbons are at match positions 1,2,3
        ("OCC(=O)CO", [1, 2, 3])  # for glycerol with one acyl substitution at one OH (we use same indices)
    ]
    backbone_found = False
    backbone_atom_indices = []  # list of tuples: (match_tuple, list_of_carbon_indices from the match)
    for pattern, carbon_positions in glycerol_smarts_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt)
        if matches:
            for m in matches:
                # m is a tuple of atom indices; we take the indices at the specified positions as the glycerol carbons.
                if max(carbon_positions) < len(m):
                    backbone_atom_indices.append( (m, [m[i] for i in carbon_positions]) )
                    backbone_found = True
            if backbone_found:
                break
    if not backbone_found:
        return False, "No glycerol backbone detected (using patterns OCC(O)CO or OCC(=O)CO)"
    
    # --- 2. Define helper function: is the oxygen part of an ester?
    def is_ester_oxygen(o_atom):
        # Check if an oxygen is attached to a carbon that is double bonded (i.e. in a carbonyl).
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon neighbor
                for bond in nbr.GetBonds():
                    if bond.GetBondTypeAsDouble() == 2.0:
                        # Found a double bond; now check that the other partner is oxygen.
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    # --- 3. Helper: calculate maximum chain length (number of carbon atoms)
    # using a DFS along C–C bonds (ignoring rings). This simple traversal will return the maximum number of connected carbons.
    def max_chain_length_from(atom, visited):
        length = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                # skip ring atoms to favor linear (aliphatic) chains:
                if nbr.IsInRing():
                    continue
                visited.add(nbr.GetIdx())
                branch_length = 1 + max_chain_length_from(nbr, visited)
                if branch_length > length:
                    length = branch_length
                visited.remove(nbr.GetIdx())
        return length

    # --- 4. Look from the glycerol backbone for an oxygen that acts as an ether linkage.
    # We require that one of the glycerol carbons (from our pattern match) is bound to an oxygen that:
    #   a) Is not part of the backbone (i.e. not in the match tuple).
    #   b) Does not look like an ester oxygen.
    #   c) Has on its "non-glycerol" side a carbon that yields a long aliphatic chain (at least 6 carbons).
    ether_found = False
    for match_tuple, glycerol_carbons in backbone_atom_indices:
        for gc_idx in glycerol_carbons:
            gc_atom = mol.GetAtomWithIdx(gc_idx)
            for nbr in gc_atom.GetNeighbors():
                # Only consider oxygen neighbors that are not part of the glycerol backbone match.
                if nbr.GetAtomicNum() != 8:
                    continue
                if nbr.GetIdx() in match_tuple:
                    continue
                # Do not consider oxygens that still carry (explicit) hydrogens (these are -OH groups)
                if nbr.GetTotalNumHs() > 0:
                    continue
                # Skip if oxygen appears to be in an ester linkage.
                if is_ester_oxygen(nbr):
                    continue
                # Now, from this candidate ether oxygen, we expect one of its neighbors (other than the glycerol carbon)
                # to be a carbon that gives rise to a long aliphatic chain.
                # (There should be exactly 2 heavy neighbors for an ether oxygen; we already used one (the glycerol carbon)).
                cell_neighbors = [nb for nb in nbr.GetNeighbors() if nb.GetIdx() != gc_idx and nb.GetAtomicNum() == 6]
                if not cell_neighbors:
                    continue
                for alkyl_c in cell_neighbors:
                    # Start a DFS from this carbon. Include this carbon as count=1.
                    chain_len = 1 + max_chain_length_from(alkyl_c, {alkyl_c.GetIdx()})
                    if chain_len >= 6:
                        ether_found = True
                        break
                if ether_found:
                    break
            if ether_found:
                break
        if ether_found:
            break

    if not ether_found:
        return False, "No non-ester (ether) linkage found attached to glycerol that leads to a long alkyl chain"
    
    return True, "Molecule contains a glycerol backbone with at least one alkyl chain attached via an ether linkage"

# Example usage: (run as a script)
if __name__ == "__main__":
    # Example test with one of the provided molecules:
    test_smiles = "CCCCCCCCCCCCCC\\C=C\\OC[C@@H](O)CO"  # 1-[(E)-hexadecen-1-yl]-sn-glycerol example
    result, reason = is_ether_lipid(test_smiles)
    print("SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)