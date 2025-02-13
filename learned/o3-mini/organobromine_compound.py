"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compound
Definition (for this classifier):
  – A compound containing at least one bromine bond that (likely) can be considered as a 
    carbon–bromine bond.
  – We accept a bond if Br is attached to an aliphatic (non‐aromatic) carbon.
  – We also accept a Br–C bond when the carbon is aromatic if that aromatic carbon appears to be
    a benzylic site (i.e. belongs to only one ring and still has at least one hydrogen).
  – We further allow Br attached to an aromatic heteroatom if that heteroatom is directly bound to at least
    one aromatic carbon with a hydrogen.
  – As a special case, we “force‐accept” N‑bromosuccinimide.
"""

from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES.
    Our heuristic acceptance criteria are:
      1. A Br attached to a carbon -> if the carbon is aliphatic, accept.
         If the carbon is aromatic, accept only if it belongs only to one ring and has at least one hydrogen.
      2. A Br attached to a non‐carbon atom is accepted only if that atom is aromatic, and it is connected
         to at least one aromatic carbon with an available hydrogen (assuming resonance delocalization).
      3. A special case: N-bromosuccinimide (NBS) is accepted.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as an organobromine compound, False otherwise.
        str: Reason for classification.
    """
    # Attempt to parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- Special case: N-bromosuccinimide ---
    # NBS is a known reagent where Br is attached to nitrogen in a cyclic imide.
    nbs_smarts = "BrN1C(=O)CCC1=O"
    nbs_pattern = Chem.MolFromSmarts(nbs_smarts)
    if mol.HasSubstructMatch(nbs_pattern):
        return True, "Matches special case: N-bromosuccinimide"

    # Iterate over all bonds to find one that might mark the compound as organobromine.
    for bond in mol.GetBonds():
        # Identify whether one end is bromine.
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetAtomicNum() == 35:
            br_atom = atom1
            partner = atom2
        elif atom2.GetAtomicNum() == 35:
            br_atom = atom2
            partner = atom1
        else:
            continue  # Neither atom is Br

        # Case 1: Br directly bound to carbon.
        if partner.GetAtomicNum() == 6:
            if not partner.GetIsAromatic():
                # Aliphatic carbon–Br; accept.
                return True, "Contains a bromine bonded to an aliphatic carbon."
            else:
                # Bromine bound to an aromatic carbon.
                # Accept only if the carbon is in a single aromatic ring and still has at least one hydrogen.
                ring_count = partner.GetRingInfo().NumRings(partner.GetIdx())
                if ring_count == 1 and partner.GetTotalNumHs() > 0:
                    return True, ("Contains a bromine bonded to an aromatic carbon "
                                  "in a benzylic-like environment (single ring, with hydrogen present).")
                else:
                    # Likely part of a fused system or heavily substituted; reject this bond.
                    continue

        # Case 2: Br bound to a non‐carbon atom.
        else:
            # First, if the partner is aromatic, we can try to see if its neighborhood suggests resonance.
            if partner.GetIsAromatic():
                # Look for at least one neighboring aromatic carbon with at least one hydrogen.
                for nbr in partner.GetNeighbors():
                    if nbr.GetIdx() == br_atom.GetIdx():
                        continue  # Skip the bromine itself.
                    if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic() and nbr.GetTotalNumHs() > 0:
                        return True, ("Bromine attached to an aromatic heteroatom that is connected "
                                      "to an aromatic carbon with available hydrogen.")
            # If partner is non-carbon and not an aromatic heteroatom (or does not meet our criteria),
            # we do not accept this bond.
            continue

    # If we iterated over all bonds without finding an acceptable pattern.
    return False, "No suitable bromine bond found."

# --- Example usage (for testing) ---
if __name__ == '__main__':
    test_smiles = [
        "[O-][N+](=O)c1ccc(CBr)cc1",                        # 4-nitrobenzyl bromide: Br at aliphatic (benzylic) C.
        "C(O)(=O)CCC/C=C(\\CC/C=C\\CCCCCCCCCCCCC(C)C)/Br",    # 6-bromo-23-methyl-tetracosa-5E,9Z-dienoic acid.
        "Brc1ccc2[nH]c3cnccc3c2c1",                          # Eudistomin N (false positive in previous attempt)
        "C1CN(CCC1(C2=CC=C(C=C2)Br)O)CCCC(=O)C3=CC=C(C=C3)F",  # Bromperidol (false positive previously)
        "BrN1C(=O)CCC1=O",                                   # N-bromosuccinimide (previously missed)
    ]
    for s in test_smiles:
        result, reason = is_organobromine_compound(s)
        print(f"SMILES: {s}\n  Classified as organobromine: {result}\n  Reason: {reason}\n")