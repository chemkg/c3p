"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: Diradylglycerol
Definition: Any lipid that is glycerol bearing two substituent groups – either acyl, alkyl, or alk-1-enyl –
at any two of the three possible positions.
Heuristic (improved):
  1. Parse the molecule from the SMILES.
  2. Look for a candidate glycerol backbone defined as a set of three carbon atoms:
       • One “central” carbon that is bonded to exactly two other carbon atoms (the two terminal carbons)
       • Each of the three carbons has at least one oxygen attached via a single bond which is not part of
         the backbone.
  3. For each oxygen directly attached to a backbone carbon, decide if it is “free” (i.e. a hydroxyl)
     or “substituted” (i.e. if it connects further to an acyl/alkyl/alk-1-enyl group). Here the
     test is simply: if the oxygen carries at least one (implicit or explicit) H then it is free.
  4. If exactly one of the three oxygens is free and the other two lack hydrogen, then return True.
     Otherwise, return False with an explanation.
"""

from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol (diacyl/di-alkyl/di-alk-1-enyl glycerol) 
    based on its SMILES string according to the heuristic:
      - The molecule must contain a glycerol backbone, that is, three core (sp3) carbons
        arranged as two terminal carbons (each bonded only to the central one) and one central carbon
        bonded to both. Each backbone carbon must carry one oxygen substituent (bonded via a single bond).
      - Exactly one of these oxygens must be free (i.e. be a hydroxyl, carrying at least one H)
        while the other two should be substituted.
        
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a diradylglycerol; False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Loop over atoms to try to find a candidate glycerol backbone.
    # We look for a "central" carbon that is connected to exactly two other carbons.
    n_atoms = mol.GetNumAtoms()
    for atom in mol.GetAtoms():
        # Check only carbon atoms
        if atom.GetSymbol() != "C":
            continue
        # Find carbon neighbors (as heavy atoms) 
        carbon_neigh = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "C"]
        if len(carbon_neigh) != 2:
            continue  # not a candidate central carbon
        central = atom
        term1, term2 = carbon_neigh

        # For a glycerol backbone, each terminal carbon should have the central as its only carbon neighbor.
        for term in (term1, term2):
            # Count carbon neighbors that are not the central.
            other_carbons = [nbr for nbr in term.GetNeighbors() if nbr.GetSymbol() == "C" and nbr.GetIdx() != central.GetIdx()]
            if len(other_carbons) != 0:
                break  # not a terminal in the backbone
        else:
            # Candidate backbone found: central and the two terminals.
            # Backbone atom indices:
            backbone_idxs = {central.GetIdx(), term1.GetIdx(), term2.GetIdx()}
            oxy_info = []  # list to store (backbone_atom_idx, oxygen_idx, is_free)
            # For each backbone carbon, search for an oxygen directly attached via a single bond
            valid_backbone = True
            for b_idx in backbone_idxs:
                carbon_atom = mol.GetAtomWithIdx(b_idx)
                # Look only at single-bonded oxygen neighbors not in the backbone.
                oxy_neighbors = []
                for bond in carbon_atom.GetBonds():
                    # We want only single bonds
                    if bond.GetBondType() != Chem.BondType.SINGLE:
                        continue
                    # Identify the neighbor atom that is not the carbon_atom.
                    nbr = bond.GetOtherAtom(carbon_atom)
                    if nbr.GetSymbol() == "O" and nbr.GetIdx() not in backbone_idxs:
                        oxy_neighbors.append((nbr, bond))
                if len(oxy_neighbors) != 1:
                    # For a proper glycerol unit we expect exactly one oxygen substituent on each backbone carbon.
                    valid_backbone = False
                    break
                nbr_atom, _ = oxy_neighbors[0]
                # Determine if the oxygen is free: i.e. it carries at least one hydrogen.
                # GetTotalNumHs returns implicit+explicit hydrogens.
                is_free = (nbr_atom.GetTotalNumHs() > 0)
                oxy_info.append((b_idx, nbr_atom.GetIdx(), is_free))
            if not valid_backbone:
                continue  # try next candidate

            # Count free vs substituted oxygens on the backbone.
            free_count = sum(1 for (_, _, is_free) in oxy_info if is_free)
            sub_count = sum(1 for (_, _, is_free) in oxy_info if not is_free)
            if free_count == 1 and sub_count == 2:
                msg = ("Found glycerol backbone (central carbon idx {}) with two substituted positions "
                       "and one free hydroxyl. Backbone oxygen details: {}").format(
                           central.GetIdx(), oxy_info)
                return True, msg
            else:
                # This candidate backbone did not show the correct pattern:
                pattern_msg = ("Candidate glycerol backbone (central atom idx {}) found, "
                               "but oxygens show freeOH: {} and substituted: {}.").format(
                                   central.GetIdx(), free_count, sub_count)
                # We continue searching for another candidate.
                # (We save pattern_msg in case no candidate is found.)
                last_reason = pattern_msg

    # If no valid backbone was found after all tests:
    return False, ("No glycerol backbone with exactly two substituted oxygen groups and one free hydroxyl "
                   "was identified. " + (last_reason if 'last_reason' in locals() else ""))

# Example usage:
if __name__ == "__main__":
    # Test with one sample SMILES (a DG structure expected to be a diradylglycerol).
    test_smiles = "C([C@@](COC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)[H])O"
    classification, reason = is_diradylglycerol(test_smiles)
    print("Classification:", classification)
    print("Reason:", reason)