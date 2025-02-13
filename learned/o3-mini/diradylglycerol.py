"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: Diradylglycerol
Definition: Any lipid that is glycerol bearing two substituent groups – either acyl, alkyl,
or alk-1-enyl – at any two of the three possible positions. In other words, a diacyl (or
di-alkyl/alk-1-enyl) glycerol should have a glycerol backbone where exactly two oxygen substituents 
are acylated/alkylated while the third hydroxyl remains free.
Heuristic (improved):
  1. Parse the molecule from the SMILES.
  2. Loop over all carbons and try to find a candidate glycerol backbone:
       • The candidate “central” carbon is sp3 and is bonded to exactly two other carbon atoms.
       • Each of those two carbons (terminals) is bonded to the central carbon as its only carbon neighbor.
  3. For each of the three backbone carbons you expect to see (at least) one oxygen attached via a single bond,
     where the oxygen is not part of the backbone. If a backbone carbon does not have any oxygen substituent,
     or has more than one, then this candidate is rejected.
  4. For each backbone oxygen substituent, decide if it is “free” (–OH) or “substituted” (esterified/etherified)
     by the following logic:
         – In many cases a free –OH will show one or more (implicit) hydrogen.
         – However, if no H is reported, we also check the connectivity: if the oxygen is attached (besides its
           backbone carbon) to another heavy atom then, if that atom is a carbon that is double-bonded to an O,
           we judge it as substituted.
         – Otherwise the oxygen is flagged as free.
  5. Finally, if exactly one of the three oxygens is considered free and the other two are substituted,
     then return True.
     
If any step fails we return False with an explanation.
"""

from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol (diacyl/di-alkyl/di-alk-1-enyl glycerol)
    based on its SMILES string according to an improved heuristic.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a diradylglycerol; False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Function to decide if an oxygen attached to a backbone carbon is free.
    def oxygen_is_free(oxy, backbone_idx):
        """
        Decide if oxygen 'oxy' is a free hydroxyl.
        First, we try the standard method: if the oxygen carries any implicit or explicit H, we call it free.
        If not, we check its other heavy-atom neighbor (if any) to see if it is part of an acyl (ester) linkage.
        """
        # Check for hydrogens: This counts both explicit and implicit H.
        if oxy.GetTotalNumHs() > 0:
            return True
        # Otherwise, get heavy neighbors other than the backbone carbon.
        heavy_neighbors = [nbr for nbr in oxy.GetNeighbors() 
                           if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != backbone_idx]
        # If no heavy neighbor beyond the backbone carbon exists, consider it free.
        if not heavy_neighbors:
            return True
        # If there is a heavy neighbor, check if it is a carbonyl carbon.
        for nbr in heavy_neighbors:
            if nbr.GetSymbol() == 'C':
                # Check if any bond from the neighbor is a double bond to an oxygen.
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetSymbol() == 'O':
                            return False  # Likely esterified (substituted)
        # Default to free if no carbonyl pattern found.
        return True

    last_reason = "No candidate glycerol backbone was found."
    # Loop over all atoms looking for a candidate for the central carbon.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        # We require the candidate central carbon to be sp3 and have exactly two carbon neighbors.
        if atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "C"]
        if len(carbon_neighbors) != 2:
            continue  # This is not a candidate central carbon.
        central = atom
        term1, term2 = carbon_neighbors

        # For the terminal carbons, require that their only carbon neighbor is the central one.
        valid_terminals = True
        for term in (term1, term2):
            # Count carbon neighbors besides the central.
            other_carbons = [nbr for nbr in term.GetNeighbors() if nbr.GetSymbol() == "C" and nbr.GetIdx() != central.GetIdx()]
            if other_carbons:
                valid_terminals = False
                break
        if not valid_terminals:
            continue

        # At this stage we have a candidate glycerol backbone: central, term1, and term2.
        backbone_idxs = {central.GetIdx(), term1.GetIdx(), term2.GetIdx()}
        oxy_info = []  # list of tuples (backbone_atom_idx, oxygen_atom_idx, is_free)
        valid_backbone = True

        # For each backbone carbon, we expect exactly one oxygen substituent (via a single bond) that is not in the backbone.
        for b_idx in backbone_idxs:
            carbon_atom = mol.GetAtomWithIdx(b_idx)
            oxygen_subs = []
            for bond in carbon_atom.GetBonds():
                # Only consider single bonds
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                nbr = bond.GetOtherAtom(carbon_atom)
                if nbr.GetSymbol() == "O" and nbr.GetIdx() not in backbone_idxs:
                    oxygen_subs.append(nbr)
            if len(oxygen_subs) != 1:
                valid_backbone = False
                last_reason = ("Backbone carbon (idx {}) does not have exactly one oxygen substituent "
                               "(found {}).".format(b_idx, len(oxygen_subs)))
                break
            oxy = oxygen_subs[0]
            # Determine if this oxygen is free (hydroxyl) or substituted (ester/ether).
            free = oxygen_is_free(oxy, b_idx)
            oxy_info.append((b_idx, oxy.GetIdx(), free))
        if not valid_backbone:
            continue

        # Count free versus substituted oxygens.
        free_count = sum(1 for (_, _, is_free) in oxy_info if is_free)
        sub_count = sum(1 for (_, _, is_free) in oxy_info if not is_free)
        if free_count == 1 and sub_count == 2:
            msg = ("Found glycerol backbone (central carbon idx {}) with two substituted positions and one free hydroxyl. "
                   "Backbone oxygen details: {}").format(central.GetIdx(), oxy_info)
            return True, msg
        else:
            last_reason = ("Candidate glycerol backbone (central carbon idx {}) found, "
                           "but oxygens show freeOH: {} and substituted: {}."
                           ).format(central.GetIdx(), free_count, sub_count)
    return False, ("No glycerol backbone with exactly two substituted oxygen groups and one free hydroxyl "
                   "was identified. " + last_reason)

# Example usage:
if __name__ == "__main__":
    # Test with one sample SMILES (a DG structure expected to be a diradylglycerol)
    test_smiles = "C([C@@](COC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)[H])O"
    classification, reason = is_diradylglycerol(test_smiles)
    print("Classification:", classification)
    print("Reason:", reason)