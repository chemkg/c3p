"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
#!/usr/bin/env python
"""
Classifies: tetradecanoate ester

A tetradecanoate ester is a fatty acid ester derived from tetradecanoic (myristic) acid.
In such a moiety the acid part is a linear, saturated chain of exactly 14 carbons (including the carbonyl carbon)
that is connected via an ester bond to an alcohol (or phenol). This implementation does not rely only on a generic
ester SMARTS but instead iterates over carbonyl centers that are candidates for ester formation. For each candidate,
it checks that the acyl chain is linear (no branching) and exactly 14 carbons long.
"""

from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule (given by its SMILES string) contains an ester 
    group whose acid fragment is exactly tetradecanoic acid (i.e. a linear, saturated 14-carbon chain).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a tetradecanoate ester is detected, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # For every carbon atom, check if it can be the carbonyl carbon of an ester.
    # In an ester, the carbonyl carbon should have:
    #   (1) a double bond to one oxygen (the carbonyl oxygen)
    #   (2) a single bond to an oxygen (the ester oxygen that connected to the alcohol part)
    #   (3) a single bond to one carbon (the start of the acyl chain)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # skip non-carbons
        # Gather information about neighbors:
        dbl_oxygens = []   # oxygen double-bonded to carbonyl carbon
        sing_oxygens = []  # oxygen singly-bonded (the ester oxygen)
        carbon_neighbors = []  # carbon substituents that are candidates for the acyl chain
        
        for bond in atom.GetBonds():
            bond_type = bond.GetBondType()
            nbr = bond.GetOtherAtom(atom)
            if bond_type == Chem.BondType.DOUBLE and nbr.GetAtomicNum() == 8:
                dbl_oxygens.append(nbr)
            elif bond_type == Chem.BondType.SINGLE:
                if nbr.GetAtomicNum() == 8:
                    sing_oxygens.append(nbr)
                elif nbr.GetAtomicNum() == 6:
                    carbon_neighbors.append(nbr)
        
        # Check that we have exactly one double-bonded oxygen and at least one ester-type oxygen.
        if len(dbl_oxygens) != 1 or len(sing_oxygens) < 1:
            continue
        
        # For a classic ester, there should be one acyl chain carbon attached via a single bond.
        # (Extra substitutions on the carbonyl carbon usually mean it is not a typical fatty acid moiety.)
        if len(carbon_neighbors) != 1:
            continue
        
        # The candidate acyl chain starts at the attached carbon.
        acyl_start = carbon_neighbors[0]
        
        # Define a helper function that follows the acyl chain along single bonds.
        # It counts the number of carbon atoms from the starting carbon until a branch is encountered.
        # The "from_atom" is the previous atom (here, the carbonyl carbon) so we do not go backwards.
        def count_linear_chain(current_atom, from_atom):
            count = 1  # count the current_atom
            prev = from_atom
            while True:
                # Get all single-bonded carbon neighbors other than the one we came from.
                nbrs = [nbr for nbr in current_atom.GetNeighbors() 
                        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev.GetIdx() and 
                           current_atom.GetBondBetweenAtoms(nbr.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
                if len(nbrs) == 1:
                    count += 1
                    prev = current_atom
                    current_atom = nbrs[0]
                else:
                    # Either no continuation or branching encountered.
                    break
            return count
        
        # The total acyl chain length = 1 (the carbonyl carbon itself) + the number of carbons in the
        # linear chain starting from acyl_start.
        total_chain_length = 1 + count_linear_chain(acyl_start, atom)
        
        # If the total chain length is exactly 14, we have a tetradecanoate unit.
        if total_chain_length == 14:
            return True, ("Contains a tetradecanoate ester moiety "
                          "(detected a linear 14-carbon acyl chain derived from myristic acid)")
    
    return False, "No tetradecanoate ester moiety detected"


# Optional test cases when running as a standalone script.
if __name__ == "__main__":
    test_smiles = [
        # These SMILES strings (from the examples) are known to include tetradecanoate ester fragments.
        ("C([C@@H](COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCC)=O)", "DG(16:0/14:0/0:0)"),
        ("CCCCCCCCCCCCC(=O)OC[C@H](O)COP(O)(=O)OC[C@@H](O)CO", "1,2-ditetradecanoyl-sn-glycerol-3-phosphate"),
        ("CCCCCCCCCCCCCCCC(=O)OCC(CO)OC(=O)CCCCCCCCCCCCCC", "Tetradecanoyl test"),
    ]
    for smi, name in test_smiles:
        valid, explanation = is_tetradecanoate_ester(smi)
        print(f"{name}: {valid} --> {explanation}")