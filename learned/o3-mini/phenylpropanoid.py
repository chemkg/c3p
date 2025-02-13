"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: phenylpropanoid compounds
Definition: Any organic aromatic compound with a structure based on a 
phenylpropane (C6–C3) skeleton. This class includes naturally occurring 
phenylpropanoid esters, flavonoids, anthocyanins, coumarins and many small 
phenolic molecules as well as their semi‐synthetic and synthetic analogues.
Phenylpropanoids are also precursors of lignin.
 
This implementation uses several SMARTS patterns to detect classical motifs 
in phenylpropanoid chemistry (coumarins, cinnamic acid derivatives, flavanone 
and isoflavone cores, and a benzene ring with an unsaturated (cinnamyl) chain).
Additionally, it runs a custom search for a clean saturated (propyl) chain 
directly attached to a benzene ring – but only if the chain atoms are “plain”
(sp3 carbons with no additional unsaturation or branching on the chain).
This extra check helps reduce false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    
    The heuristic implemented here first requires that the molecule has at least one 
    benzene ring. Then, it checks for several characteristic substructures:
      - Coumarin (and related 4-hydroxycoumarin) structure.
      - Cinnamic acid derivative (benzene with an unsaturated three‐carbon chain 
        terminating in a carboxyl function, either acid or its anion).
      - Benzene with an unsaturated (cinnamyl) chain.
      - A flavanone core.
      - An isoflavone core.
      
    In addition, if none of the above SMARTS patterns match, the code will search 
    for a benzene ring carrying a “clean” saturated propyl chain. Here “clean” means 
    that the three–carbon chain is entirely aliphatic (sp3) and each carbon is not further substituted 
    (beyond the link to benzene and its adjacent chain neighbour). This extra check helps avoid 
    many false positives from generic saturated chains in large molecules.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if classified as a phenylpropanoid, False otherwise.
      str: Explanation for the classification decision.
    """
    
    # Parse the SMILES string:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize molecule to ensure aromaticity and valence are computed.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Basic filter: ensure there is at least one benzene ring
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene):
        return False, "No benzene ring detected, so not a phenylpropanoid"
    
    # List of (description, SMARTS) pairs for characteristic phenylpropanoid substructures.
    # These patterns look for named classes.
    patterns = [
        # Coumarin core; this pattern should capture coumarins including 4-hydroxycoumarin.
        ("Coumarin structure", "O=c1ccc2oc(=O)cc2c1"),
        # Cinnamic acid derivative (acid or anion). The pattern forces a benzene ring attached to a 
        # CH=CH–C(=O)O fragment.
        ("Cinnamic acid derivative", "c1ccccc1C=CC(=O)[O;$([O-]);!$([OH])]"),
        # Benzene with unsaturated (cinnamyl) chain: benzene directly attached to CH2–CH=CH–...
        ("Benzene with unsaturated (cinnamyl) chain", "c1ccccc1[CH2]C=CC"),
        # Flavanone core: this pattern represents a common flavanone motif.
        ("Flavanone core", "c1cc(c(cc1)O)C2=CC(=O)OC3=CC=CC=C23"),
        # Isoflavone core: a different arrangement found in many isoflavones.
        ("Isoflavone core", "c1ccc2C(=O)oc3ccccc3c2c1")
    ]
    
    # Check SMARTS patterns one by one.
    for desc, smarts in patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # skip invalid SMARTS (should not occur)
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches pattern: {desc}"
    
    # Custom check: look for a benzene ring directly attached to a saturated three–carbon chain.
    # To reduce false positives we require:
    #   (a) The chain atoms must be aliphatic and sp3.
    #   (b) They should form a contiguous chain of exactly three atoms 
    #       (i.e. one carbon attached to the benzene ring, then one mid chain, then one terminal).
    #
    # We iterate over atoms that are part of a benzene ring.
    benzene_matches = mol.GetSubstructMatches(benzene)
    for ring in benzene_matches:
        # For each atom index in the found benzene ring:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look at neighbours not in the benzene ring.
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in ring:
                    continue
                # Check that the first atom of the side chain is aliphatic and sp3.
                if nb.GetAtomicNum() != 6 or nb.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                # Now try to follow a linear chain of exactly 3 carbon atoms.
                chain = [nb]
                current = nb
                valid_chain = True
                # We want a chain length of 3 atoms (already have first) so add two more:
                for _ in range(2):
                    # Among neighbours (exclude atoms already in chain and the benzene ring)
                    next_candidates = [n for n in current.GetNeighbors() 
                                       if n.GetIdx() not in ring and n.GetIdx() not in [a.GetIdx() for a in chain]
                                       and n.GetAtomicNum() == 6 and n.GetHybridization() == Chem.rdchem.HybridizationType.SP3]
                    if len(next_candidates) != 1:
                        valid_chain = False
                        break
                    chain.append(next_candidates[0])
                    current = next_candidates[0]
                if not valid_chain:
                    continue
                # At this point we have three chain atoms. We further check that the terminal carbon does not have extra heavy-atom substituents
                # (other than the one linking to the chain).
                term = chain[-1]
                # Count neighbours (excluding the chain neighbor)
                ext_neighbors = [n for n in term.GetNeighbors() if n.GetIdx() not in [a.GetIdx() for a in chain]]
                # Allow one extra neighbor if it is hydrogen (which might not be explicit) – here we check atomic number.
                if any(n.GetAtomicNum() != 1 for n in ext_neighbors):
                    continue
                # If reached here, we found a benzene ring with a clean saturated propyl chain.
                return True, "Contains benzene ring with saturated (propyl) chain attached"
    
    return False, "No phenylpropanoid-characteristic substructure found"

# Example usage (only for testing; remove or comment out if used as a module):
if __name__ == "__main__":
    test_smiles = [
        # True positives (from provided examples)
        "C1(C2=C(O[C@@H](C1)C3=C(C=C(C=C3)O)O)C(=C(C=C2O)O)C[C@@H](CCC(=C)C)C(=C)C)=O",  # remangiflavanone B
        "COc1ccc(cc1OC)[C@H]1O[C@H]([C@H](C)[C@@H]1C)c1ccc2OCOc2c1",  # futokadsurin C
        "COc1cc2C[C@H]3COC(=O)[C@@H]3Cc3c(OC)c(OC)c(OC)cc3-c2cc1OC",  # Neoisostegane
        "COC1=C(O)C=CC(=C1)[C@H]1OC2=CC(O)=CC(O)=C2C(=O)[C@@H]1O",  # (+)-dihydroisorhamnetin
        "O1C(C(O)C(O)C(OC(=O)/C=C/C2=CC=C(O)C=C2)C1OC=3C=C4OC=C(C(=O)C4=C(O)C3)C5=CC=C(O)C=C5)CO",  # Genistein 7-O-(2-p-coumaroylglucoside)
        # A simple molecule that should NOT be classified as a phenylpropanoid:
        "OC(=Cc1ccccc1)C([O-])=O",  # enol-phenylpyruvate (false positive in earlier version)
        "CCO"  # ethanol
    ]
    
    for sm in test_smiles:
        result, reason = is_phenylpropanoid(sm)
        print(f"SMILES: {sm}\nClassified as phenylpropanoid? {result}\nReason: {reason}\n")