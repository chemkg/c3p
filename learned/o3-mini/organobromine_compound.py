"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compound
Definition:
  – A compound containing at least one carbon-bromine bond.
  – Accept if Br is attached to an aliphatic carbon.
  – If Br is attached to an aromatic carbon, accept only if that carbon is in exactly one aromatic ring and has at least one hydrogen.
  – Accept also if Br is attached to an aromatic heteroatom that is directly bound to at least one aromatic carbon with a hydrogen.
  – Special case: N‑bromosuccinimide is accepted.
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES.
    Our heuristic acceptance criteria are:
      1. A bromine attached to a carbon:
         - If the carbon is aliphatic, accept.
         - If the carbon is aromatic, accept only if it is part of exactly one ring and has at least one hydrogen.
      2. A bromine attached to a non‐carbon atom:
         - Accept if that atom is aromatic and has at least one neighboring aromatic carbon with an available hydrogen.
      3. Special case: N‑bromosuccinimide is accepted.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as an organobromine compound, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Special case: N-bromosuccinimide.
    # SMARTS for NBS: bromine bound to a nitrogen in a cyclic imide.
    nbs_smarts = "BrN1C(=O)CCC1=O"
    nbs_pattern = Chem.MolFromSmarts(nbs_smarts)
    if mol.HasSubstructMatch(nbs_pattern):
        return True, "Matches special case: N-bromosuccinimide"

    # Get ring info once from the molecule.
    ring_info = mol.GetRingInfo()

    # Iterate over all bonds to look for a bromine bond.
    for bond in mol.GetBonds():
        # Check if one of the atoms is bromine (atomic number 35).
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetAtomicNum() == 35:
            br_atom = atom1
            partner = atom2
        elif atom2.GetAtomicNum() == 35:
            br_atom = atom2
            partner = atom1
        else:
            continue  # Neither atom is bromine
        
        # Case 1: Bromine directly attached to a carbon.
        if partner.GetAtomicNum() == 6:
            if not partner.GetIsAromatic():
                # If the carbon is not aromatic (aliphatic), accept.
                return True, "Contains bromine bonded to an aliphatic carbon."
            else:
                # If the carbon is aromatic, check that it is part of exactly one ring and has at least one hydrogen.
                ring_count = ring_info.NumAtomRings(partner.GetIdx())
                # GetTotalNumHs() gives the total count of hydrogens (implicit and explicit).
                if ring_count == 1 and partner.GetTotalNumHs() > 0:
                    return True, ("Contains bromine bonded to an aromatic carbon in a benzylic-like environment "
                                  "(single ring with available hydrogen).")
                else:
                    # Not acceptable due to being in multiple rings or lacking hydrogens.
                    continue
        # Case 2: Bromine attached to a non–carbon atom.
        else:
            # Accept only if the partner atom is aromatic.
            if partner.GetIsAromatic():
                # Check if at least one neighboring aromatic carbon (exclude the bromine neighbor) has at least one hydrogen.
                for nbr in partner.GetNeighbors():
                    if nbr.GetIdx() == br_atom.GetIdx():
                        continue
                    if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic() and nbr.GetTotalNumHs() > 0:
                        return True, ("Bromine attached to an aromatic heteroatom that is connected to an aromatic carbon "
                                      "with available hydrogen.")
            # If not matching our criteria, continue to other bonds.
            continue
    
    # After checking all bonds, if no suitable bond was found, reject.
    return False, "No suitable bromine bond found."

# --- Example usage (for testing) ---
if __name__ == '__main__':
    test_smiles = [
        "[O-][N+](=O)c1ccc(CBr)cc1",                        # 4-nitrobenzyl bromide
        "C(O)(=O)CCC/C=C(\\CC/C=C\\CCCCCCCCCCCCC(C)C)/Br",    # 6-bromo-23-methyl-tetracosa-5E,9Z-dienoic acid
        "[C@@]123O[C@]([C@H](C)CC1(C)C)(CC(O[C@@]([C@@H](C)O)(CC(OC(C2)[C@H](C)[C@](O3)([C@H](CC[C@@H](C=4C(=CC(=C(C4)O)Br)Br)OC)C)[H])=O)[H])=O)O",  # 19-bromoaplysiatoxin
        "BrN1C(=O)CCC1=O",                                   # N-bromosuccinimide special case
    ]
    
    for s in test_smiles:
        result, reason = is_organobromine_compound(s)
        print(f"SMILES: {s}\n  Classified as organobromine: {result}\n  Reason: {reason}\n")