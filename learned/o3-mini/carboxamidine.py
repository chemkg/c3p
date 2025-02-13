"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: Carboxamidine containing compounds
Definition: Compounds having the structure RC(=NR)NR2 (commonly the -C(=NH)NH2 group, or variants thereof).
This improved version uses a broader SMARTS pattern and inspects the immediate bonding environment 
to avoid flagging peptide-like groups.
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group based on its SMILES string.
    The group is defined as RC(=NR)NR2 (for example, the common -C(=NH)NH2 group) and related variants.
    
    The algorithm is:
      1. Parse the SMILES.
      2. Look for a carboxamidine fragment using a more inclusive SMARTS:
         “[CX3](=[NX2,NX3,NX4+])[NX3,NX4]”
         This should catch neutral and charged variants.
      3. For each SMARTS match, inspect the environment of the central carbon.
         The carbon should be connected to exactly two nitrogens (part of the amidine) and one extra substituent.
         If that extra substituent is directly a “peptide‐like” fragment (e.g. a carbon with a double bond to oxygen),
         then the match is flagged as part of a peptide environment.
      4. Additionally, if the overall molecule contains peptide‐bond fragments (using “N-C(=O)C”) and is large,
         we flag the molecule as peptide‐like.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a carboxamidine group is detected in a valid (nonpeptide-like) context, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a broad SMARTS pattern to capture carboxamidine groups.
    # This pattern accepts variations in the bonding (including the possibility of formal charges)
    # It requires a carbon (3-connected) with a double bond to a nitrogen (which may have 2,3 or even 4 connections)
    # and a single bond to a second nitrogen.
    pattern_amidine = Chem.MolFromSmarts("[CX3](=[NX2,NX3,NX4+])[NX3,NX4]")
    matches = mol.GetSubstructMatches(pattern_amidine)
    if not matches:
        return False, "No carboxamidine group detected."
    
    # Define a simple SMARTS for a peptide bond fragment.
    peptide_pattern = Chem.MolFromSmarts("N-C(=O)C")
    # If the molecule is large and contains peptide bonds, flag it as peptide-like.
    is_peptide_mol = mol.HasSubstructMatch(peptide_pattern) and mol.GetNumHeavyAtoms() > 30
    
    valid_match_found = False
    
    # Iterate over all matches for the carboxamidine group:
    for match in matches:
        # By our SMARTS the first atom (match[0]) is the central carbon.
        c_idx = match[0]
        c_atom = mol.GetAtomWithIdx(c_idx)
        
        # We expect the central carbon to have exactly three neighbors.
        # Two of these should be the amidine nitrogens (contained in the match).
        neighbor_indices = [nbr.GetIdx() for nbr in c_atom.GetNeighbors()]
        extra_neighbors = [idx for idx in neighbor_indices if idx not in match[1:]]
        
        peptide_like_env = False
        # If there is an extra substituent, inspect it:
        for idx in extra_neighbors:
            nbr_atom = mol.GetAtomWithIdx(idx)
            # Check each bond from the substituent:
            for bond in nbr_atom.GetBonds():
                # If the bond is a double bond and the bonded atom is oxygen (carbonyl-like),
                # mark as peptide-like.
                if bond.GetBondTypeAsDouble() == 2.0:
                    other = bond.GetOtherAtom(nbr_atom)
                    if other.GetAtomicNum() == 8:
                        peptide_like_env = True
                        break
            if peptide_like_env:
                break
        
        # Accept the match if no peptide-like extra substituent is found.
        if not peptide_like_env:
            valid_match_found = True
            break

    if valid_match_found and not is_peptide_mol:
        return True, "Carboxamidine group detected via substructure match."
    else:
        if is_peptide_mol:
            return False, "WRONGLY CLASSIFIED Carboxamidine group detected via substructure match (peptide-like environment)."
        else:
            return False, "No valid carboxamidine group detected after environmental inspection."

# Example usage (for testing):
if __name__ == "__main__":
    test_cases = {
        # True positives (from provided examples)
        "4-[5-(3,5-dichlorophenyl)-5-(trifluoromethyl)-4,5-dihydro-1,2-oxazol-3-yl]-N-[(methoxyamino)methylidene]-2-methylbenzamide":
            "CON=CNC(=O)c1ccc(cc1C)C1=NOC(C1)(c1cc(Cl)cc(Cl)c1)C(F)(F)F",
        "acetamiprid":
            "CN(Cc1ccc(Cl)nc1)C(C)=NC#N",
        "A-234 nerve agent":
            "P(=O)(N=C(N(CC)CC)C)(OCC)F",
        "tegaserod":
            "CCCCCNC(=N)N\\N=C\\c1c[nH]c2ccc(OC)cc12",
        "formamidine":
            "[H]C(N)=N",
        "benzamidine":
            "NC(=N)c1ccccc1",
        "cyflufenamid":
            "Fc1ccc(c(\\C(NC(=O)Cc2ccccc2)=N\\OCC2CC2)c1F)C(F)(F)F",
        "2-formamido-N(1)-(5-phospho-D-ribosyl)acetamidine":
            "O[C@H]1[C@@H](O)C(NC(=N)CNC=O)O[C@@H]1COP(O)(O)=O",
        "pyrantel":
            "CN1CCCN=C1\\C=C\\c1cccs1",
        "debrisoquin":
            "NC(=N)N1CCc2ccccc2C1",
        "amitraz":
            "CN(C=Nc1ccc(C)cc1C)C=Nc1ccc(C)cc1C",
        "N'-(3-hydroxyphenyl)-N,N-dimethylformamidine":
            "CN(C)C=Nc1cccc(O)c1",
        "A-232 nerve agent":
            "P(=O)(N=C(N(CC)CC)C)(OC)F",
        "4-hydroxydebrisoquin":
            "NC(=N)N1CC(O)c2ccccc2C1",
        "2,5-bis{[4-(N-ethylamidino)]phenyl}furan":
            "CC\\N=C(/N)c1ccc(cc1)-c1ccc(o1)-c1ccc(cc1)C(\\N)=N\\CC",
        "Acrylamidine":
            "N=C(N)C=C",
        "8-(pyrimidin-2-ylamino)naphthalene-2-carboximidamide":
            "NC(=N)c1ccc2cccc(Nc3ncccn3)c2c1",
        "dabigatran etexilate":
            "CCCCCCOC(=O)\\N=C(\\N)c1ccc(NCc2nc3cc(ccc3n2C)C(=O)N(CCC(=O)OCC)c2ccccn2)cc1",
        "5-hydroxyectoine":
            "CC1=N[C@@H]([C@@H](O)CN1)C(O)=O",
        "N'-(2,4-dimethylphenyl)-N-methylformamidine":
            "CN=CNc1ccc(C)cc1C",
        "guanidine":
            "NC(N)=N",
        "lofexidine":
            "CC(Oc1c(Cl)cccc1Cl)C1=NCCN1",
        "7-formamidino-7-deazaguanine":
            "NC(=N)c1c[nH]c2nc(N)[nH]c(=O)c12",
        "acetamidine":
            "CC(N)=N",
        "(Z)-acetamiprid":
            "CN(Cc1ccc(Cl)nc1)C(\\C)=N/C#N",
        # False positives (peptide-like environments)
        "Arg-Ser-Tyr":
            "O=C(N[C@H](C(=O)N[C@@H](CC1=CC=C(O)C=C1)C(O)=O)CO)[C@@H](N)CCCN=C(N)N",
        "N(omega)-methyl-L-argininate":
            "CNC(=N)NCCC[C@H](N)C([O-])=O",
        "4-guanidinobutanoic acid":
            "C(=N)(N)NCCCC(O)=O",
        "L-Acetopine":
            "OC(=O)C(NCC(O)=O)CCCN=C(N)N",
        # False negatives (should be carboxamidine but missed previously)
        "melagatran":
            "[H][C@]1(CCN1C(=O)[C@H](NCC(O)=O)C1CCCCC1)C(=O)NCc1ccc(cc1)C(N)=N",
        "nizatidine":
            "CNC(NCCSCc1csc(CN(C)C)n1)=C[N+]([O-])=O",
        "kasugamycin":
            "C[C@H]1O[C@H](O[C@@H]2[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]2O)[C@@H](N)C[C@@H]1NC(=N)C(O)=O"
    }
    
    for name, smi in test_cases.items():
        res, explanation = is_carboxamidine(smi)
        print(f"Test: {name}\nSMILES: {smi}\nResult: {res}\nReason: {explanation}\n{'-'*60}")