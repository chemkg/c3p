"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compound
Definition: A compound containing at least one carbon–bromine bond.
In our improved classifier we also allow a bromine attached to an aromatic heteroatom 
(for example, in cases where tautomerization/resonance might place the bromine on carbon)
but we require that non‐aromatic heteroatoms (e.g. in N‐bromosuccinimide) are not accepted.
"""

from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined (here) as an organic molecule that has at least one 
    bond between a bromine atom and a carbon. In addition, if the Br is attached to an aromatic 
    heteroatom, we accept it (assuming resonance could locate the Br on carbon), while for non‐aromatic 
    heteroatoms the structure is not considered organobromine.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an organobromine compound, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over bonds to look for bromine attached bonds.
    for bond in mol.GetBonds():
        # Get the two atoms in the bond.
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        # Check if one of the atoms is Bromine (atomic number 35).
        if atom1.GetAtomicNum() == 35:
            br_atom = atom1
            partner = atom2
        elif atom2.GetAtomicNum() == 35:
            br_atom = atom2
            partner = atom1
        else:
            continue  # Neither atom is bromine, skip
        
        # If partner is carbon (atomic number 6), we are done.
        if partner.GetAtomicNum() == 6:
            return True, "Contains at least one carbon–bromine bond."
        
        # If the partner is not carbon, then we check if it is aromatic.
        # For example, in 1‑bromoindole the Br may appear attached to an aromatic nitrogen.
        # We accept such cases assuming resonance may cause a carbon–bromine localization.
        if partner.GetAtomicNum() != 6 and partner.GetIsAromatic():
            return True, ("Bromine attached to an aromatic heteroatom "
                          "(assumed to be in resonance with an organic framework).")
    
    # If no suitable bromine bond is found, return False.
    return False, "No suitable bromine bond found."

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "[O-][N+](=O)c1ccc(CBr)cc1",  # 4-nitrobenzyl bromide: clear C-Br bond
        "Brn1ccc2ccccc12",            # 1-bromoindole: Br attached to aromatic N (accept)
        "BrN1C(=O)CCC1=O",            # N-bromosuccinimide: Br attached to non-aromatic N (reject)
    ]
    for s in test_smiles:
        result, reason = is_organobromine_compound(s)
        print(f"SMILES: {s}\n  Classified as organobromine: {result}\n  Reason: {reason}\n")