"""
Classifies: CHEBI:76578 diradylglycerol
"""
#!/usr/bin/env python3
"""
Classifies: Diradylglycerol
Definition: Any lipid that is glycerol bearing two substituent groups – either acyl, alkyl,
or alk-1-enyl – at any two of the three possible positions. That is, a diacyl (or di-alkyl/alk-1-enyl)
glycerol should have a glycerol backbone in which exactly two of its oxygen groups are esterified (or etherified)
while the third remains a free hydroxyl.
Heuristic:
  1. Parse the molecule from the SMILES string.
  2. Loop over all atoms to find a candidate glycerol backbone.
     For a candidate central carbon:
       - It must be sp3 and have exactly two carbon neighbors.
       - Each of its two attached (terminal) carbons must have the central carbon as its only carbon neighbor.
  3. For each of the three glycerol backbone carbons (central + two terminals),
     look for exactly one oxygen substituent (via a single bond) that is not itself part of the backbone.
  4. For each oxygen in the backbone substituents, classify it as “free” (–OH) or “substituted” (acyl/alkyl)
     via the following improved logic:
         • If it has any attached hydrogen, call it free.
         • Otherwise, check the non-backbone neighbors.
             – If any neighbor is phosphorus, reject this candidate (phospholipid-like).
             – If a neighbor is carbon that bears a double bond to another oxygen (a carbonyl) then mark as substituted.
         • If no hydrogen is found and no convincing carbonyl is found, assume it is substituted.
  5. Finally, if exactly one oxygen is free and two are substituted, we classify the input as a diradylglycerol.
If no candidate passes, return False with a summary reason.
"""

from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol (a diacyl/di-alkyl/di-alk-1-enyl glycerol)
    based on its SMILES string using a heuristic for a glycerol backbone and analysis of oxygen substituents.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a diradylglycerol; False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved function to decide if an oxygen atom (attached to a backbone carbon) is free (OH) or substituted.
    def oxygen_is_free(oxy, backbone_idx):
        """
        Determine if the oxygen 'oxy' is a free hydroxyl (returns True) or is substituted (returns False).
        Additionally, if the oxygen is attached to a phosphorus atom, we mark this candidate as not a DG.
        """
        # Use explicit and implicit hydrogen count.
        num_H = oxy.GetNumExplicitHs() + oxy.GetNumImplicitHs()
        if num_H > 0:
            # Likely a free hydroxyl.
            return True
        
        # Examine non-backbone neighbors.
        for nbr in oxy.GetNeighbors():
            if nbr.GetIdx() == backbone_idx:
                continue
            # If the oxygen is attached to phosphorus, then it belongs to e.g. a phosphate headgroup.
            if nbr.GetAtomicNum() == 15:
                # Signal an ambiguous situation by returning None.
                return None
            if nbr.GetSymbol() == 'C':
                # Look at bonds from the neighbor to see if it bears a carbonyl.
                for bond in nbr.GetBonds():
                    # Identify a double bond to another oxygen from this neighbor.
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetSymbol() == 'O' and other.GetIdx() != oxy.GetIdx():
                            return False  # This oxygen is part of an ester linkage.
        # No hydrogen and no clear carbonyl detected: assume substituted (e.g. an ether or vinyl linkage).
        return False
    
    last_reason = "No candidate glycerol backbone was found."
    # Loop over all atoms to search for a candidate central carbon of a glycerol backbone.
    for atom in mol.GetAtoms():
        # We require the candidate central atom to be carbon.
        if atom.GetSymbol() != "C":
            continue
        # The candidate must be sp3 hybridized.
        if atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        # The candidate central carbon should have exactly two carbon neighbors.
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "C"]
        if len(carbon_neighbors) != 2:
            continue
        central = atom
        term1, term2 = carbon_neighbors

        # For each terminal carbon, its only carbon neighbor should be the central carbon.
        valid_terminals = True
        for term in (term1, term2):
            others = [nbr for nbr in term.GetNeighbors() if nbr.GetSymbol() == "C" and nbr.GetIdx() != central.GetIdx()]
            if others:
                valid_terminals = False
                break
        if not valid_terminals:
            continue

        # We have a candidate glycerol backbone (central, term1, term2).
        backbone_idxs = {central.GetIdx(), term1.GetIdx(), term2.GetIdx()}
        oxy_info = []  # Each entry: (backbone_atom_idx, oxygen_atom_idx, is_free)
        valid_backbone = True

        # For each backbone carbon we expect exactly one oxygen connected by a single (non-backbone) bond.
        for b_idx in backbone_idxs:
            carbon_atom = mol.GetAtomWithIdx(b_idx)
            oxygen_subs = []
            for bond in carbon_atom.GetBonds():
                # Only consider single bonds.
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                nbr = bond.GetOtherAtom(carbon_atom)
                # Exclude atoms that belong to the backbone.
                if nbr.GetSymbol() == "O" and nbr.GetIdx() not in backbone_idxs:
                    oxygen_subs.append(nbr)
            if len(oxygen_subs) != 1:
                valid_backbone = False
                last_reason = ("Backbone carbon (idx {}) does not have exactly one oxygen substituent (found {})."
                               .format(b_idx, len(oxygen_subs)))
                break
            oxy = oxygen_subs[0]
            free = oxygen_is_free(oxy, b_idx)
            # If our oxygen assessment returns None, then the candidate backbone is part of a phospholipid.
            if free is None:
                valid_backbone = False
                last_reason = ("Backbone oxygen (idx {}) attached to carbon (idx {}) is connected to phosphorus."
                               .format(oxy.GetIdx(), b_idx))
                break
            oxy_info.append((b_idx, oxy.GetIdx(), free))
        if not valid_backbone:
            continue

        # In a diradylglycerol we expect exactly one free hydroxyl and two substituted oxygens.
        free_count = sum(1 for (_, _, is_free) in oxy_info if is_free)
        sub_count  = sum(1 for (_, _, is_free) in oxy_info if not is_free)
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