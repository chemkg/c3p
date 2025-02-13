"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: ether lipid
Definition: A lipid similar in structure to a glycerolipid but in which one or more of the carbon atoms
on glycerol is bonded to an alkyl chain via an ether (C–O–C) linkage instead of an ester (CO–O) linkage.
Our strategy is to first detect a glycerol backbone (using SMARTS for free glycerol or monoacylated glycerol)
but then require that the backbone carbons are not part of any ring. Then for each backbone carbon we search
its oxygen neighbors. We require that an oxygen candidate (a) is not part of the backbone, (b) does not have
any explicit hydrogens (so that it is an ether substituent, not a hydroxyl) and (c) is not part of an ester linkage.
Finally the candidate oxygen must lead into a long aliphatic chain (at least 6 connected carbons via non‐ring C–C bonds).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    
    The algorithm:
      1. Parse and sanitize the SMILES.
      2. Look for a glycerol-like backbone using two SMARTS patterns: 
             a) free glycerol "OCC(O)CO" 
             b) monoacylated glycerol "OCC(=O)CO".
         In each match we require that the three matched backbone carbons (positions 1,2,3 in the pattern)
         are not in a ring.
      3. For each backbone carbon, check its oxygen neighbors. For each oxygen we:
             - Ensure it is not part of the backbone match.
             - Require that it carries no explicit hydrogen (so that it is already substituted, indicating an ether bond).
             - Check that it is not an ester oxygen (i.e. not bound to a carbonyl carbon).
      4. For a candidate oxygen, examine its non-backbone neighbor (which should be a carbon) and use a DFS
         (while ignoring ring atoms) to count the length of the attached carbon chain.
         If at least 6 carbon atoms are connected linearly, we count it as the ether alkyl chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an ether lipid, False otherwise.
        str: Explanation message.
    """
    # --- 1. Parse and sanitize the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Sanitization failed: " + str(e)
    
    # --- 2. Search for a glycerol-like backbone.
    # We define two SMARTS patterns: one for free glycerol (HO–CH2–CHOH–CH2–OH) and
    # one for monoacylated glycerol (HO–CH2–CHOH–CH2–O–C(=O)-). Using the pattern "OCC(O)CO" and "OCC(=O)CO".
    # For both, we assume that the backbone carbons come from the 2nd, 3rd and 4th atoms in the match.
    glycerol_smarts_patterns = [
        ("OCC(O)CO", [1, 2, 3]),   # free glycerol backbone; ex: HO–CH2–CHOH–CH2–OH
        ("OCC(=O)CO", [1, 2, 3])   # monoacylated backbone (one acyl substitution)
    ]
    backbone_found = False
    backbone_matches = []  # list of tuples: (match_tuple, list of backbone carbon indices)
    for pattern, carbon_positions in glycerol_smarts_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt)
        if matches:
            for m in matches:
                # Ensure our match is long enough and that the identified backbone carbons are not in any ring.
                if max(carbon_positions) >= len(m):
                    continue
                backbone_atoms = [m[i] for i in carbon_positions]
                if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in backbone_atoms):
                    # skip matches in which any backbone atom is in a ring (to avoid sugars and cyclic polyols)
                    continue
                backbone_matches.append((m, backbone_atoms))
                backbone_found = True
            if backbone_found:
                break
    if not backbone_found:
        return False, "No appropriate glycerol backbone detected (using patterns OCC(O)CO or OCC(=O)CO, acyclic match)"
    
    # --- 3. Helper function: Check if an oxygen is part of an ester (i.e. attached to a carbonyl).
    def is_ester_oxygen(o_atom):
        # For a candidate oxygen, look at its carbon neighbors.
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon neighbor
                for bond in nbr.GetBonds():
                    # check for a double bond from the carbon to oxygen (i.e. a carbonyl)
                    if bond.GetBondTypeAsDouble() == 2.0:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    # --- 4. Helper function: DFS to count maximum chain length from a carbon atom along C–C bonds.
    def max_chain_length_from(atom, visited):
        max_length = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                # Ignore ring atoms to favor linear unbranched chains.
                if nbr.IsInRing():
                    continue
                visited.add(nbr.GetIdx())
                branch_length = 1 + max_chain_length_from(nbr, visited)
                if branch_length > max_length:
                    max_length = branch_length
                visited.remove(nbr.GetIdx())
        return max_length

    # --- 5. Look from any glycerol backbone match for an oxygen that fulfills our ether criteria.
    ether_found = False
    for match_tuple, glycerol_carbons in backbone_matches:
        # For each glycerol carbon in this backbone, look at neighbors.
        for gc_idx in glycerol_carbons:
            gc_atom = mol.GetAtomWithIdx(gc_idx)
            for nbr in gc_atom.GetNeighbors():
                # Only consider oxygen neighbors.
                if nbr.GetAtomicNum() != 8:
                    continue
                # Skip oxygens that are part of the glycerol backbone.
                if nbr.GetIdx() in match_tuple:
                    continue
                # We require that the oxygen is not a free hydroxyl;
                # if it still carries explicit hydrogens, then it is not an ether substitution.
                if nbr.GetTotalNumHs() > 0:
                    continue
                # Also skip if this oxygen is in a ring: many sugars have cyclic ethers.
                if nbr.IsInRing():
                    continue
                # Do not consider oxygens in an ester linkage.
                if is_ester_oxygen(nbr):
                    continue
                # Now check: the oxygen should be bound to a carbon (other than the glycerol carbon) 
                # that leads to a long aliphatic chain.
                candidate_carbons = [n for n in nbr.GetNeighbors() if n.GetIdx() != gc_idx and n.GetAtomicNum() == 6]
                if not candidate_carbons:
                    continue
                for alkyl_c in candidate_carbons:
                    # Start a DFS from this carbon. Count this carbon as 1.
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
        return False, "No non-ester (ether) linkage found attached to an acyclic glycerol backbone that leads to a long alkyl chain"
    
    return True, "Molecule contains an acyclic glycerol backbone with at least one alkyl chain attached via an ether (C–O–C) linkage"

# Example usage: (run as a script)
if __name__ == "__main__":
    # To test, try one of the provided examples. Here we test with 1-[(E)-hexadecen-1-yl]-sn-glycerol.
    test_smiles = "CCCCCCCCCCCCCC\\C=C\\OC[C@@H](O)CO"
    result, reason = is_ether_lipid(test_smiles)
    print("SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)