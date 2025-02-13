"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: Carboxamidine containing compounds
Definition: Compounds having the structure RC(=NR)NR2 (commonly the -C(=NH)NH2 group,
or variants thereof).  This version uses a refined SMARTS pattern and inspects the bonding
environment to reduce false positives from peptide-like moieties.
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group (RC(=NR)NR2) in a valid,
    non-peptide-like context.

    The algorithm is:
      1. Parse the SMILES string.
      2. Identify possible amidine groups using a broad SMARTS pattern that captures variations,
         including charged species.
      3. For every match, check:
           a. The central carbon (first atom in the match) must have exactly three bonds.
           b. Two of these bonds should be to the amidine nitrogens (from the SMARTS match),
              and the remaining one is the substituent (R).
           c. If the extra substituent is a carbon that is directly double-bonded to an oxygen,
              then it is likely in a peptide-like environment and the match is rejected.
      4. In addition, if the molecule overall contains typical peptide bond fragments
         (e.g. "N-C(=O)C") and has a large skeleton, the molecule is considered peptide-like.
      5. If at least one match passes and the overall molecule does not appear peptide-like,
         we classify the molecule as containing a carboxamidine group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if a valid carboxamidine group is detected (and not in a peptide-like environment),
              False otherwise.
        str: Explanation of the decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a relatively broad SMARTS pattern for amidine groups.
    # This pattern matches a central sp2 carbon with a double bond to a nitrogen (which can be neutral or charged)
    # and a single bond to a second nitrogen.
    pattern_amidine = Chem.MolFromSmarts("[CX3](=[NX2,NX3,NX4+])[NX3,NX4]")
    matches = mol.GetSubstructMatches(pattern_amidine)
    if not matches:
        return False, "No carboxamidine group detected."
    
    # Define a simple SMARTS for a peptide bond fragment.
    peptide_pattern = Chem.MolFromSmarts("N-C(=O)C")
    # If the molecule contains peptide bonds and is relatively large, flag it as peptide-like.
    is_peptide_mol = mol.HasSubstructMatch(peptide_pattern) and mol.GetNumHeavyAtoms() > 30
    
    valid_match_found = False
    # Iterate over each match for the amidine group.
    for match in matches:
        # Expect the match to have exactly 3 atoms: [central C, N1, N2]
        if len(match) != 3:
            continue
        c_idx = match[0]
        c_atom = mol.GetAtomWithIdx(c_idx)
        
        # The central carbon should be three-coordinate (attached to exactly 3 atoms).
        if c_atom.GetDegree() != 3:
            continue
        
        # Get all neighboring atoms for the central carbon.
        neighbor_indices = [nbr.GetIdx() for nbr in c_atom.GetNeighbors()]
        # Identify the extra substituent that is not a nitrogen in the amidine group.
        extra_neighbors = [idx for idx in neighbor_indices if idx not in match[1:]]
        
        # There should be exactly one such substituent.
        if len(extra_neighbors) != 1:
            continue
        
        extra_idx = extra_neighbors[0]
        extra_atom = mol.GetAtomWithIdx(extra_idx)
        
        # Check if the extra substituent is a carbon that is directly double-bonded to oxygen.
        peptide_like_env = False
        if extra_atom.GetAtomicNum() == 6:
            # Look at all bonds emanating from the substituent atom.
            for bond in extra_atom.GetBonds():
                # If the bond is a double bond (C=O) and the other atom is oxygen, mark as peptide-like.
                if bond.GetBondTypeAsDouble() == 2.0:
                    other = bond.GetOtherAtom(extra_atom)
                    if other.GetAtomicNum() == 8:
                        peptide_like_env = True
                        break

        if peptide_like_env:
            continue  # Skip this match as it appears to be in a peptide-like context.
        
        # If we get here, the match appears to be a valid carboxamidine fragment.
        valid_match_found = True
        break

    if valid_match_found and not is_peptide_mol:
        return True, "Carboxamidine group detected via substructure match."
    else:
        if is_peptide_mol:
            return False, "WRONGLY CLASSIFIED Carboxamidine group detected via substructure match (peptide-like environment)."
        else:
            return False, "No valid carboxamidine group detected after environmental inspection."

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_cases = {
        # True positives (examples from provided list)
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
        print(f"Test: {name}")
        print(f"SMILES: {smi}")
        print(f"Result: {res}")
        print(f"Reason: {explanation}")
        print("-" * 60)