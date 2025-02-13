"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: Monoradylglycerol
Definition:
   Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position.
   
Our approach:
 1. Parse the SMILES and add hydrogens.
 2. Look for any three-carbon contiguous chain (sp3 carbons with the central atom connected only to the two others in the chain)
    as a candidate glycerol backbone.
 3. For each candidate backbone, collect oxygen substituents (neighbors outside the backbone) on each carbon.
    In “free” glycerol, each carbon would have a free hydroxyl group (an oxygen with at least one hydrogen).
 4. In a monoradylglycerol exactly two of these oxygens should be free and one oxygen should be substituted (no hydrogens).
 5. For the substituted oxygen, follow its bond away from the backbone and check that a carbon chain of length ≥ 6 is attached.
 6. If any candidate passes these criteria, return True with an explanation.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines whether a molecule is a monoradylglycerol.
    A monoradylglycerol is defined as any glycerol (a three‐carbon backbone with three oxygen substituents)
    in which exactly one of the three oxygens is substituted by a fatty (acyl, alkyl, or alk-1-enyl) chain,
    while the remaining two are free –OH groups.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a monoradylglycerol, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that free -OH groups carry hydrogen numbers.
    mol = Chem.AddHs(mol)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"
    
    # Helper: recursively compute maximum carbon chain length from a given carbon.
    def longest_chain_length(atom, visited):
        if atom.GetAtomicNum() != 6:
            return 0
        max_length = 1
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                length = 1 + longest_chain_length(nbr, visited.copy())
                if length > max_length:
                    max_length = length
        return max_length

    # Find candidate glycerol backbones.
    # We require: three sp3 carbons connected in a chain;
    # the middle carbon must connect to both terminal carbons, and the terminal ones should be linked only to the middle (within the chain).
    sp3_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3]
    visited_chains = set()
    
    for a in sp3_carbons:
        for b in a.GetNeighbors():
            if b.GetAtomicNum() != 6 or b.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            for c in b.GetNeighbors():
                if c.GetAtomicNum() != 6 or c.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                if c.GetIdx() == a.GetIdx():
                    continue  # avoid loops
                # Order the candidate backbone by indices to avoid repeats.
                backbone = tuple(sorted([a.GetIdx(), b.GetIdx(), c.GetIdx()]))
                if backbone in visited_chains:
                    continue
                visited_chains.add(backbone)
                # Check connectivity: b should be connected to both a and c.
                if not(a in b.GetNeighbors() and c in b.GetNeighbors()):
                    continue
                # Check that a and c are terminal in the candidate (i.e. their only neighbor among backbone atoms is b)
                if sum(1 for nbr in a.GetNeighbors() if nbr.GetIdx() in backbone) != 1:
                    continue
                if sum(1 for nbr in c.GetNeighbors() if nbr.GetIdx() in backbone) != 1:
                    continue
                
                # Now, for each of these three backbone carbons, collect oxygen substituents (neighbors not in backbone).
                backbone_atoms = [a, b, c]
                oxygens_per_carbon = []
                valid_candidate = True
                for carbon in backbone_atoms:
                    oxygens = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [at.GetIdx() for at in backbone_atoms]]
                    if len(oxygens) != 1:
                        valid_candidate = False
                        break
                    oxygens_per_carbon.append((carbon, oxygens[0]))
                if not valid_candidate:
                    continue
                
                # Now count free –OH versus substituted oxygen.
                free_count = 0
                substituted_info = []  # (parent_carbon, oxygen atom)
                for parent, oxy in oxygens_per_carbon:
                    # A free hydroxyl should have at least one hydrogen after adding explicit hydrogens.
                    if oxy.GetTotalNumHs() > 0:
                        free_count += 1
                    else:
                        substituted_info.append((parent, oxy))
                # For a monoradylglycerol exactly 2 free hydroxyls and 1 substituted oxygen are expected.
                if free_count != 2 or len(substituted_info) != 1:
                    continue

                # Check the fatty (substituted) portion.
                parent, sub_oxy = substituted_info[0]
                # Find neighbors of the substituted oxygen that are not part of the glycerol backbone.
                fatty_candidates = [nbr for nbr in sub_oxy.GetNeighbors() if nbr.GetIdx() not in [at.GetIdx() for at in backbone_atoms]]
                if not fatty_candidates:
                    continue
                fatty_found = False
                for nbr in fatty_candidates:
                    # If the neighbor is a carbon, it might be directly attached (alkyl/alk-1-enyl)
                    # or it might be a carbonyl carbon (acyl chain).
                    if nbr.GetAtomicNum() == 6:
                        # Check if this carbon is a carbonyl (has a double-bonded oxygen neighbor).
                        bonds = [mol.GetBondBetweenAtoms(nbr.GetIdx(), nb.GetIdx()) for nb in nbr.GetNeighbors() if nb.GetAtomicNum() == 8 and nb.GetIdx() != sub_oxy.GetIdx()]
                        is_carbonyl = any(bond.GetBondTypeAsDouble() == 2 for bond in bonds if bond is not None)
                        if is_carbonyl:
                            # In an acyl chain, follow the chain from the carbonyl carbon: look for a carbon neighbor (not sub_oxy)
                            for acyl_nbr in nbr.GetNeighbors():
                                if acyl_nbr.GetIdx() == sub_oxy.GetIdx():
                                    continue
                                if acyl_nbr.GetAtomicNum() == 6:
                                    if longest_chain_length(acyl_nbr, set()) >= 6:
                                        fatty_found = True
                                        break
                        else:
                            # For a typical alkyl/alkenyl chain.
                            if longest_chain_length(nbr, set()) >= 6:
                                fatty_found = True
                    if fatty_found:
                        break
                if not fatty_found:
                    continue
                
                # If we reach here, we have a candidate glycerol backbone with 2 free -OH and one fatty substituent.
                return True, "Found glycerol backbone with exactly one fatty substituent (monoradylglycerol)"
    
    return False, "Could not identify a glycerol backbone with a single fatty substituent"

# Example usage and testing:
if __name__ == "__main__":
    test_examples = [
        # These examples should be monoradylglycerols.
        "O(C(=O)CCCCCCCCC/C=C/CCCCCC)CC(O)CO",     # MG(18:1(11E)/0:0/0:0)[rac]
        "O(C(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)C(CO)CO",  # MG(0:0/22:5(...)/0:0)
        "C(CCCCCCCCCCCCC(=O)OCC(O)CO)",            # simplified 1-monooctanoylglycerol-like structure
        "CCCC(=O)OCC(O)CO",                        # monobutyrin (should fail fatty chain check)
    ]
    for s in test_examples:
        result, reason = is_monoradylglycerol(s)
        print("SMILES:", s)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 50)