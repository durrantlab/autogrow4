import json
from typing import List, Optional, Set, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem
import re
import gzip
import random
from rdkit.Chem.Draw import rdMolDraw2D
import os
import sys
import enum
from rdkit.Chem import Draw

class ReactionStatus(enum.Enum):
    """Enum to represent the outcome of a reaction test."""
    SUCCESS = "SUCCESS"
    FORWARD_FAILURE = "FORWARD_FAILURE"
    REVERSE_FAILURE = "REVERSE_FAILURE"
    INVALID_REACTANTS = "INVALID_REACTANTS"  # New status for reactants that don't match group SMARTS
    REVERSE_SKIPPED = "REVERSE_SKIPPED"  # New status for skipped reverse reactions due to exceptions

def check_reverse_reaction_exception(reactants: List[str], exception_smarts: str, reaction_name: str, reactant_idx: int) -> bool:
    """
    Check if any of the reactants matches the reverse reaction exception SMARTS pattern.
    Returns True if an exception is found (should skip reverse reaction), False otherwise.
    """
    try:
        # Convert exception SMARTS to molecule pattern
        exception_mol = Chem.MolFromSmarts(exception_smarts)
        if exception_mol is None:
            print(f"WARNING: Invalid reverse_reaction_exception SMARTS pattern: {exception_smarts} in {reaction_name}, set {reactant_idx}")
            return False
        
        # Check each reactant against the exception pattern
        for i, reactant_smiles in enumerate(reactants):
            reactant_mol = Chem.MolFromSmiles(reactant_smiles)
            if reactant_mol is None:
                continue  # Skip invalid reactants (they'll be caught elsewhere)
            
            if reactant_mol.HasSubstructMatch(exception_mol):
                print(f"REVERSE SKIP: Reactant {i+1} ({reactant_smiles}) matches reverse_reaction_exception pattern ({exception_smarts}) in {reaction_name}, set {reactant_idx}")
                return True
        
        return False
        
    except Exception as e:
        print(f"ERROR: Failed to process reverse_reaction_exception for {reaction_name}, set {reactant_idx}: {e}")
        return False

# --- NEW FUNCTION TO VALIDATE REACTANTS AGAINST GROUP SMARTS ---
def validate_reactants_against_smarts(reactants: List[str], group_smarts: List[str], reaction_name: str, reactant_idx: int) -> bool:
    """
    Check if the reactants match their corresponding group SMARTS patterns.
    Returns True if all reactants match, False otherwise.
    """
    if len(reactants) != len(group_smarts):
        print(f"WARNING: Reactant count ({len(reactants)}) doesn't match group SMARTS count ({len(group_smarts)}) for {reaction_name}, set {reactant_idx}")
        return False
    
    for i, (reactant_smiles, smarts_pattern) in enumerate(zip(reactants, group_smarts)):
        # Convert reactant SMILES to molecule
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
        if reactant_mol is None:
            print(f"INVALID: Could not create molecule from SMILES: {reactant_smiles} in {reaction_name}, set {reactant_idx}")
            return False
        
        # Convert SMARTS pattern to molecule pattern
        smarts_mol = Chem.MolFromSmarts(smarts_pattern)
        if smarts_mol is None:
            print(f"INVALID: Could not create pattern from SMARTS: {smarts_pattern} in {reaction_name}, set {reactant_idx}")
            return False
        
        # Check if the reactant matches the pattern
        if not reactant_mol.HasSubstructMatch(smarts_mol):
            print(f"SKIP: Reactant {i+1} ({reactant_smiles}) does not match group SMARTS pattern ({smarts_pattern}) in {reaction_name}, set {reactant_idx}")
            return False
    
    return True

# --- NEW FUNCTION TO GENERATE DETAILED DEBUGGING PROMPT ---
def generate_llm_debug_prompt(
    reaction_name: str,
    form_name: str,
    forward_smarts: str,
    forward_reactants: List[str],
    forward_products: List[List[str]],
    reverse_smarts_list: List[str],
    reverse_products_raw: Set[str],
    original_reactants_cleaned: Set[str],
    reverse_products_cleaned: Set[str],
):
    """
    Formats a detailed, copy-paste-friendly debug prompt for an LLM
    when a reverse reaction fails.
    """
    missing_reactants = original_reactants_cleaned - reverse_products_cleaned
    
    prompt = f"""
################################################################################
# LLM DEBUGGING PROMPT: REVERSE REACTION FAILED
################################################################################

Hello, I'm debugging a chemical reaction definition and need your help. A reverse reaction failed to regenerate the original reactants. Here is a complete summary of the failure:

---
1. REACTION CONTEXT
---
- Reaction Name: "{reaction_name}"
- Reactant Form Tested: {form_name}

---
2. FORWARD REACTION (This part appears to be working)
---
- Forward Reaction SMARTS:
`{forward_smarts}`

- These Reactants were provided:
{forward_reactants}

- The Forward Reaction produced this set of products:
{forward_products}

---
3. REVERSE REACTION (This is where the failure occurred)
---
- The products from the forward reaction were then fed into the following Reverse Reaction SMARTS:
{reverse_smarts_list}

- The complete, raw set of molecules generated by the reverse reactions was:
{reverse_products_raw}

---
4. THE FAILURE
---
- My script expected the original reactants to be a subset of the molecules produced by the reverse reaction.
- Cleaned Original Reactants (Expected):
{original_reactants_cleaned}

- Cleaned Reverse Products (Actual):
{reverse_products_cleaned}

- CRITICAL: The following molecules were expected but were NOT found in the reverse reaction products:
{missing_reactants}

---
5. QUESTION
---
Based on the reactants, the forward products, and the failed reverse products, can you identify the error in the reverse reaction SMARTS list and provide a corrected version that would successfully regenerate the missing reactants?

################################################################################
"""
    print(prompt)


def validate_functional_group_defs(reactions):
    func_grp_defs = {}
    for reaction_name, reaction_info in reactions.items():
        for i, func_grp_name in enumerate(reaction_info["functional_groups"]):
            func_grp_smarts = reaction_info["group_smarts"][i]
            if func_grp_name not in func_grp_defs:
                func_grp_defs[func_grp_name] = func_grp_smarts
            else:
                # Check if the existing definition matches the new one
                if func_grp_defs[func_grp_name] != func_grp_smarts:
                    print(f"Warning: Functional group '{func_grp_name}' has multiple definitions:")
                    print(f"  Existing: {func_grp_defs[func_grp_name]}")
                    print(f"  New:      {func_grp_smarts}")
                    # Optionally, you could raise an error here instead of just printing
                    raise ValueError(f"Multiple definitions for functional group '{func_grp_name}'")

# --- NEW VALIDATION FUNCTION TO CHECK SMARTS CONSISTENCY ---
def validate_smarts_consistency_and_examples(reactions: dict) -> None:
    """
    Validates that 'group_smarts' and 'reaction_string' are semantically
    equivalent and that the example reactants match these definitions.
    """
    print("\n--- Running SMARTS Consistency and Example Validation ---")
    for reaction_name, reaction_info in reactions.items():
        if "SKIP" in reaction_name:
            continue
        
        num_reactants = reaction_info.get("num_reactants")
        if not num_reactants:
            continue

        reaction_string = reaction_info["reaction_string"]
        reactant_side = reaction_string.split(">>")[0]
        reaction_reactant_smarts_list = reactant_side.split('.')

        # --- Basic Sanity Checks ---
        if len(reaction_reactant_smarts_list) != num_reactants:
            raise ValueError(
                f"FATAL: In '{reaction_name}', 'reaction_string' implies "
                f"{len(reaction_reactant_smarts_list)} reactants, but 'num_reactants' is {num_reactants}."
            )
        if len(reaction_info["group_smarts"]) != num_reactants:
             raise ValueError(
                f"FATAL: In '{reaction_name}', 'group_smarts' has "
                f"{len(reaction_info['group_smarts'])} entries, but 'num_reactants' is {num_reactants}."
            )

        for i in range(num_reactants):
            group_smarts = reaction_info["group_smarts"][i]
            reaction_reactant_smarts = reaction_reactant_smarts_list[i]
            example_reactant_smiles = reaction_info["example_rxn_reactants"][i]
            
            # 1. Clean the reaction SMARTS by removing atom mapping
            reaction_reactant_smarts_cleaned = re.sub(r':\d+', '', reaction_reactant_smarts)

            # 2. Perform semantic comparison using RDKit
            mol_from_group = Chem.MolFromSmarts(group_smarts)
            mol_from_reaction = Chem.MolFromSmarts(reaction_reactant_smarts_cleaned)

            if mol_from_group is None:
                raise ValueError(f"FATAL: Invalid SMARTS in 'group_smarts' for '{reaction_name}', reactant {i+1}: {group_smarts}")
            if mol_from_reaction is None:
                raise ValueError(f"FATAL: Invalid SMARTS in 'reaction_string' for '{reaction_name}', reactant {i+1}: {reaction_reactant_smarts_cleaned}")

            canonical_group = Chem.MolToSmarts(mol_from_group)
            canonical_reaction = Chem.MolToSmarts(mol_from_reaction)
            
            if canonical_group != canonical_reaction:
                error_msg = f"""
FATAL: SMARTS definition mismatch in reaction '{reaction_name}' for reactant {i+1}.
  - The 'group_smarts' definition is:
    '{group_smarts}' (Canonical: '{canonical_group}')
  - The reactant in 'reaction_string' is:
    '{reaction_reactant_smarts_cleaned}' (Canonical: '{canonical_reaction}')
These should be semantically equivalent but are not. Please correct the definitions in reactions.json.
"""
                raise ValueError(error_msg)

            # 3. Validate that the example reactant matches both patterns
            example_mol = Chem.MolFromSmiles(example_reactant_smiles)
            if example_mol is None:
                 raise ValueError(f"FATAL: Invalid SMILES for example reactant {i+1} in '{reaction_name}': {example_reactant_smiles}")

            # if not example_mol.HasSubstructMatch(mol_from_group):
            #     raise ValueError(
            #         f"FATAL: Example reactant validation failed for reaction '{reaction_name}', reactant {i+1}.\n"
            #         f"  - Example SMILES: '{example_reactant_smiles}'\n"
            #         f"  - This molecule does not contain the required pattern from 'group_smarts':\n"
            #         f"    '{group_smarts}'\n"
            #         f"Please provide a correct example reactant."
            #     )
            # We only need to check one, since we proved they are equivalent above
            # if not example_mol.HasSubstructMatch(mol_from_reaction): ...

    print("--- Validation Successful: All definitions are consistent. ---\n")


def clean_undefined_stereo(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # Check if there are any undefined stereocenters (bonds with ~)
    if '~' in smiles:
        # Remove only undefined stereochemistry, keep defined stereo
        Chem.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol)
    else:
        # No undefined stereo, return as-is
        return smiles
    
def force_bond_resolution(smiles):
    if "~" not in smiles:
        return smiles
    
    # First try: just remove the ~ symbol (most common solution)
    try:
        cleaned_smiles = smiles.replace("~", "")
        mol = Chem.MolFromSmiles(cleaned_smiles)
        if mol is not None:
            Chem.SanitizeMol(mol)
            # Let RDKit figure out aromaticity
            Chem.SetAromaticity(mol)
            return Chem.MolToSmiles(mol)
    except Exception as e:
        pass
    
    # Fallback: replace with single bond
    try:
        cleaned_smiles = smiles.replace("~", "-")
        mol = Chem.MolFromSmiles(cleaned_smiles)
        if mol is not None:
            Chem.SanitizeMol(mol)
            return Chem.MolToSmiles(mol)
    except Exception as e:
        pass
    
    return None

def clean_up_smiles(smiles: str) -> Optional[str]:
    """A simple function to clean up common verbose SMILES representations."""
    # TODO: Do these need to be applied to SMILES in the source databases?
    smiles = smiles.replace("[CH2]", "C")
    smiles = smiles.replace("[CH3]", "C")
    smiles = smiles.replace("[CH]", "C")
    smiles = smiles.replace("[CH2-]", "C")
    smiles = smiles.replace("[CH-]", "C")
    smiles = smiles.replace("[C]", "C")
    smiles = smiles.replace("[c]", "c")
    smiles = smiles.replace("[OH-]", "O")
    smiles = smiles.replace("[OH+]", "O")
    smiles = smiles.replace("[IH]", "I")  # Why this?
    smiles = smiles.replace("[NH2+]", "N")
    smiles = smiles.replace("[N]", "N")
    smiles = smiles.replace("[NH2]", "N")
    smiles = smiles.replace("[NH]", "N")
    smiles = smiles.replace("[nH+][n-]", "[nH]n")  # This is a common pattern in azides
    smiles = smiles.replace("N#[N+][N-]", "[N-]=[N+]=N")
    smiles = smiles.replace("[N-][N+]#N", "N=[N+]=[N-]")
    smiles = smiles.replace("C[N+]=C=O", "CN=C=O")
    smiles = smiles.replace("[B]", "B")  # Remove brackets around boron
    # smiles = smiles.replace("[O-]", "O")
    smiles = smiles.replace("[SH-]", "S")
    smiles = smiles.replace("[Cl+]", "Cl")
    smiles = smiles.replace("[Br+]", "Br")
    smiles = smiles.replace("[I+]", "I")
    smiles = smiles.replace("(=[O+])", "(=O)")  # Remove + from carbonyls

    # Strange tetrazole resonance structure I see in some SMILES
    smiles = smiles.replace("C2=N[N+]=NN2", "C2=NN=NN2")
    # Add more replacements if needed

    # Also, smiles should be canonicalized
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # print("Invalid SMILES: ", smiles) # Silenced
        return None
    
    try:
        # Let RDKit handle hydrogen assignment properly
        Chem.SanitizeMol(mol)
        # Remove explicit hydrogens and let RDKit recalculate
        mol = Chem.RemoveHs(mol)
        smiles = Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        # print(f"Sanitization failed for {smiles}: {e}") # Silenced
        return None
    
    # smiles = Chem.MolToSmiles(mol, canonical=True)

    smiles = clean_undefined_stereo(smiles)
    smiles = force_bond_resolution(smiles)

    # if "~" in smiles:
        # import pdb; pdb.set_trace()  # Debugging point to check if undefined stereo is still present

    # Note that this makes invalid SMILES sometimes. But necessary to identify
    # equivalent structures in some cases.
    # smiles = smiles.replace("[nH]", "n")  

    return smiles

# --- IMPORTANT NOTE ON REACTION PRODUCT CLEANUP ---
# It is an inevitable and expected outcome that applying general reverse reaction SMARTS
# to a complex product will generate some chemically nonsensical structures. This happens
# because SMARTS is a pattern-matching language, not a chemistry engine. It will
# apply the transformation rule wherever the pattern matches, even if it results in
# an invalid molecule (e.g., breaking a bond inside a stable aromatic ring).
#
# When these malformed structures are passed to RDKit for sanitization (`clean_up_smiles`),
# RDKit correctly identifies them as invalid and fails to produce a valid SMILES string.
# This is not a bug in the code or the reaction definitions, but rather a feature of
# this robust validation approach. The framework correctly filters out these "garbage"
# products and only considers the valid ones for the reverse reaction check.
#
# Therefore, failures to clean up or sanitize certain products are silently ignored.
# Do not attempt to "fix" the reverse reaction strings to prevent this, as it would
# require making them overly specific and brittle, which is a worse problem. This
# behavior should be accepted as part of the process.

def run_rxn(rxn_string: str, 
            reactants_smi_list: List[str], 
            reaction_name: str, 
            reactant_idx: int, 
            form_name: str,
            reaction_step: str) -> List[List[str]]:
    """
    Runs a reaction and returns only valid, sanitized product sets.
    Failures to sanitize certain products are silently ignored, as they are
    expected artifacts of applying general SMARTS rules.
    """
    rxn = AllChem.ReactionFromSmarts(rxn_string)
    
    # Convert all reactant SMILES strings to Mol objects
    reactants_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smi_list]
    
    # Check if any molecule conversion failed
    if any(mol is None for mol in reactants_mols):
        # Handle error, e.g., by logging a warning and returning empty
        print(f"Warning: Could not create one or more reactant molecules from SMILES: {reactants_smi_list}")
        return []

    # Run the reaction. RunReactants expects a sequence (tuple or list) of Mol objects.
    product_sets = rxn.RunReactants(reactants_mols)
    
    predicted_products = []
    for product_set in product_sets:
        product_smiles = []
        for product in product_set:
            if product is not None:
                try:
                    # Attempt to get a SMILES string first
                    raw_smiles = Chem.MolToSmiles(product)
                    # Now attempt to clean it
                    cleaned_smiles = clean_up_smiles(raw_smiles)

                    if cleaned_smiles:
                        product_smiles.append(cleaned_smiles)
                    # Silently ignore products that fail to clean up. This is expected.
                
                except Exception as e:
                    # RDKit can fail even on MolToSmiles for some malformed products.
                    # Silently ignore these as they are expected artifacts.
                    pass

        if product_smiles:
            # Sort the smiles in each set to make comparisons order-independent
            predicted_products.append(sorted(product_smiles))
            
    return predicted_products

# --- NEWLY ADDED VALIDATION FUNCTION ---
def validate_reaction_examples(reactions: dict) -> None:
    """
    Validates that the 'example_rxn_reactants' for each reaction, when processed
    by 'reaction_string', produce the 'example_rxn_product'.
    """
    print("\n--- Running Reaction Example Validation ---")
    for reaction_name, reaction_info in reactions.items():
        if "SKIP" in reaction_name:
            continue

        # Get the required fields
        example_reactants = reaction_info.get("example_rxn_reactants")
        expected_product_smi = reaction_info.get("example_rxn_product")
        reaction_smarts = reaction_info.get("reaction_string")

        # Skip if the example is not fully defined
        if not all([example_reactants, expected_product_smi, reaction_smarts]):
            print(f"INFO: Skipping example validation for '{reaction_name}' due to missing 'example_rxn_reactants', 'example_rxn_product', or 'reaction_string'.")
            continue

        # Clean the expected product for a fair comparison
        cleaned_expected_product = clean_up_smiles(expected_product_smi)
        if not cleaned_expected_product:
            raise ValueError(
                f"FATAL: Could not process 'example_rxn_product' for reaction '{reaction_name}': "
                f"'{expected_product_smi}'"
            )

        # Apply the same cleaning to the reactants
        cleaned_reactants = [clean_up_smiles(r) for r in example_reactants]
        if None in cleaned_reactants:
             raise ValueError(
                f"FATAL: Could not process one or more 'example_rxn_reactants' for reaction '{reaction_name}': "
                f"{example_reactants}"
            )

        # Run the forward reaction using the provided example
        actual_product_sets = run_rxn(
            reaction_smarts,
            cleaned_reactants,
            reaction_name,
            -1,  # Using -1 to indicate this is a validation run
            "example_validation",
            "forward"
        )

        # Flatten the list of product sets into a single set of unique products
        actual_products_flat = {product for product_set in actual_product_sets for product in product_set}

        # Check if the cleaned expected product is in the set of actual products
        if cleaned_expected_product not in actual_products_flat:
            error_msg = f"""
FATAL: Example reaction validation failed for '{reaction_name}'.
  - Reactants Used: {cleaned_reactants}
  - Reaction SMARTS: '{reaction_smarts}'
  - Expected Product (Canonical): '{cleaned_expected_product}'
  - Actual Products Generated: {actual_products_flat or 'None'}
The expected product was NOT found in the generated products. Please check the reaction definition and example in reactions.json.
"""
            raise ValueError(error_msg)

    print("--- Validation Successful: All reaction examples work as expected. ---\n")


# --- REVISED apply_prereaction FUNCTION ---
def apply_prereaction(smiles: str,
                      prereaction_strings: List[str],
                      reaction_name: str,
                      reactant_idx: int) -> str:
    """
    Apply a list of prereaction strings to standardize a molecule.
    If any prereaction yields products, return the first product from the first successful reaction.
    If no prereaction yields products, return the original SMILES.
    
    Args:
        smiles: The input SMILES string
        prereaction_strings: List of SMARTS reaction strings to try
        reaction_name: Name of the reaction for context
        reactant_idx: Index of the reactant set for context
        
    Returns:
        Standardized SMILES string
    """
    for prereaction_string in prereaction_strings:
        try:
            # Try to run this prereaction on the molecule, providing all required context
            products = run_rxn(
                prereaction_string, [smiles],
                reaction_name, reactant_idx, "prereaction", "prereaction"
            )
            
            # If we got products, use the first product from the first product set
            if products and len(products) > 0 and len(products[0]) > 0:
                standardized_smiles = products[0][0]  # First product from first set
                print(f"  Applied prereaction: {smiles} -> {standardized_smiles}")
                return standardized_smiles
                
        except Exception as e:
            # If this prereaction fails, try the next one
            print(f"  Prereaction failed: {prereaction_string} on {smiles}: {e}")
            continue
    
    # If no prereaction worked, return original
    # print(f"  No prereaction applied to: {smiles}")
    return smiles

# --- REVISED get_forward_reactant_groups FUNCTION ---
def get_forward_reactant_groups(functional_groups: List[str],
                                reaction_name: str,
                                debug_prereactions: Optional[List[str]] = None) -> List[List[str]]:
    reactant_groups = []

    for functional_group in functional_groups:
        flnm = f"./groups/{functional_group}.smi.gz"
        try:
            with gzip.open(flnm, "rt") as f:
                lines = f.readlines()

                # Pick 100 random reactants from the group
                random.shuffle(lines)
                lines = lines[:100]

                lines = [Chem.MolToSmiles(Chem.MolFromSmiles(line.strip().split("\t")[0])) for line in lines if line.strip()]

                reactant_groups.append(lines)
        except FileNotFoundError:
            print(f"Warning: Functional group file {flnm} not found. Skipping this group.")
            continue

    reactant_groups2 = []
    if not reactant_groups or not any(reactant_groups):
        return []
        
    min_len = min(len(g) for g in reactant_groups)
    for i in range(min_len):
        # First, do the forward reaction
        forward_reactants = [clean_up_smiles(reactant_groups[j][i]) for j in range(len(reactant_groups))]
        forward_reactants = [f for f in forward_reactants if f]
        
        # Apply debug prereactions if they exist
        if debug_prereactions:
            print(f"Applying debug prereactions to reactant set {i}:")
            standardized_reactants = []
            for j, reactant in enumerate(forward_reactants):
                # Apply prereactions to each reactant individually
                standardized_reactant = apply_prereaction(
                    reactant, debug_prereactions, reaction_name, i
                )
                standardized_reactants.append(standardized_reactant)
            forward_reactants = standardized_reactants
        
        reactant_groups2.append(forward_reactants)

    return reactant_groups2

def save_reaction_visualization(reaction_name: str, forward_reactants: List[str], forward_products: List[List[str]], 
                               all_reverse_products: List[str], original_reactants_set: set, 
                               reverse_products_set: set, reactant_index: int, 
                               form_name: str, success: bool = True) -> None:
    """
    Save a visualization of the PASS/FAIL reaction pathway.
    """
    viz_dir = "reaction_visualizations"
    os.makedirs(viz_dir, exist_ok=True)

    # Create a safe filename
    safe_reaction_name = re.sub(r'[^\w\-_\.]', '_', reaction_name)
    
    # Get actual molecules for original reactants
    reactant_mols = [Chem.MolFromSmiles(smi) for smi in forward_reactants if Chem.MolFromSmiles(smi) is not None]
    product_mols = [Chem.MolFromSmiles(smi) for p_set in forward_products for smi in p_set if Chem.MolFromSmiles(smi)]
    reverse_mols = [Chem.MolFromSmiles(smi) for smi in reverse_products_set if Chem.MolFromSmiles(smi)]

    all_mols = reactant_mols + [None] + product_mols + [None] + reverse_mols
    legends = ([f"Reactant {i+1}" for i in range(len(reactant_mols))] + ["→"] + 
              [f"Fwd Prod {i+1}" for i in range(len(product_mols))] + ["Rev→"] + 
              [f"Rev Prod {i+1}" for i in range(len(reverse_mols))])

    status = "PASS" if success else "FAIL"
    png_filename = f"{viz_dir}/{safe_reaction_name}_{reactant_index:03d}_{form_name}_{status}.png"
    
    try:
        max_per_row = max(len(reactant_mols), len(product_mols), len(reverse_mols)) + 1
        img = Draw.MolsToGridImage(all_mols, molsPerRow=max_per_row, subImgSize=(200, 200), legends=legends, useSVG=False)
        img.save(png_filename)
    except Exception as e:
        print(f"Failed to save PNG visualization {png_filename}: {e}")

# --- MODIFIED execute_and_verify_reaction FUNCTION ---
def execute_and_verify_reaction(
    reaction_name: str,
    reaction_info: dict,
    reactants_to_test: List[str],
    reactant_idx: int,
    form_name: str,
) -> ReactionStatus:
    """
    Runs the full forward and reverse reaction cycle and returns a status.
    Now includes validation that reactants match their group SMARTS patterns
    and checks for reverse reaction exceptions.
    """
    # NEW: Validate reactants against group SMARTS before proceeding
    group_smarts = reaction_info.get("group_smarts", [])
    if not validate_reactants_against_smarts(reactants_to_test, group_smarts, reaction_name, reactant_idx):
        return ReactionStatus.INVALID_REACTANTS
    
    # NEW: Check for reverse reaction exception
    reverse_exception_smarts = reaction_info.get("reverse_reaction_exception")
    if reverse_exception_smarts:
        if check_reverse_reaction_exception(reactants_to_test, reverse_exception_smarts, reaction_name, reactant_idx):
            return ReactionStatus.REVERSE_SKIPPED
    
    reaction_string = reaction_info["reaction_string"]
    reverse_reaction_strings = reaction_info.get("reverse_reaction_strings", [])

    # Forward reaction call now includes context for visualization
    forward_products = run_rxn(
        reaction_string, reactants_to_test, 
        reaction_name, reactant_idx, form_name, "forward"
    )

    if not forward_products:
        print(f"INFO: No valid products generated for forward reaction ({form_name}): {reactants_to_test}")
        return ReactionStatus.FORWARD_FAILURE

    # Reverse reaction
    all_reverse_products = []
    for reverse_reactants_set in forward_products:
        for reverse_reaction_string in reverse_reaction_strings:
            # Reverse reaction call also includes context
            reverse_products = run_rxn(
                reverse_reaction_string, reverse_reactants_set,
                reaction_name, reactant_idx, form_name, "reverse"
            )
            for r in reverse_products:
                all_reverse_products.extend(iter(r))
    
    # Comparison logic remains the same
    original_reactants_set = {clean_up_smiles(p) for p in reactants_to_test if clean_up_smiles(p)}
    reverse_products_set = {clean_up_smiles(p) for p in all_reverse_products if clean_up_smiles(p)}

    if original_reactants_set.issubset(reverse_products_set):
        print(f"Passed! Reaction: {reaction_name}, Form: {form_name}")
        save_reaction_visualization(reaction_name, reactants_to_test, forward_products,
                                  all_reverse_products, original_reactants_set,
                                  reverse_products_set, reactant_idx, form_name, success=True)
        return ReactionStatus.SUCCESS
    else:
        # Save visualization for a true reverse failure
        save_reaction_visualization(reaction_name, reactants_to_test, forward_products,
                                  all_reverse_products, original_reactants_set,
                                  reverse_products_set, reactant_idx, form_name, success=False)
        # And print the LLM debug prompt
        generate_llm_debug_prompt(
            reaction_name=reaction_name,
            form_name=form_name,
            forward_smarts=reaction_string,
            forward_reactants=reactants_to_test,
            forward_products=forward_products,
            reverse_smarts_list=reverse_reaction_strings,
            reverse_products_raw=set(all_reverse_products),
            original_reactants_cleaned=original_reactants_set,
            reverse_products_cleaned=reverse_products_set,
        )
        return ReactionStatus.REVERSE_FAILURE

# --- MAIN EXECUTION LOGIC ---

with open("reactions.json", "r") as file:
    contents = file.read()
    # Replace some placeholders
    contents = contents.replace("AZIDE_SEARCH", "[$(N=[N+]=[N-]),$([N-][N+]#N)]")
    reactions = json.loads(contents)

# --- RUN ALL VALIDATION CHECKS BEFORE TESTING ---
validate_functional_group_defs(reactions)
# validate_smarts_consistency_and_examples(reactions)
validate_reaction_examples(reactions) # <-- New validation function is called here

# --- GENERATE FUNCTIONAL GROUPS OUTPUT FILE ---
def generate_functional_groups_output(reactions: dict) -> None:
    """
    Generate a functional_groups_out.json file containing all unique functional group definitions.
    """
    functional_groups_dict = {}
    
    for reaction_name, reaction_info in reactions.items():
        if "SKIP" in reaction_name:
            continue
            
        functional_groups = reaction_info.get("functional_groups", [])
        group_smarts = reaction_info.get("group_smarts", [])
        
        # Pair up functional group names with their SMARTS definitions
        for func_group, smarts in zip(functional_groups, group_smarts):
            if func_group not in functional_groups_dict:
                functional_groups_dict[func_group] = smarts
            else:
                # Verify consistency (should already be validated above)
                if functional_groups_dict[func_group] != smarts:
                    print(f"WARNING: Inconsistent SMARTS for {func_group}")
    
    # Sort by functional group name for consistent output
    sorted_groups = dict(sorted(functional_groups_dict.items()))
    
    # Write to file with nice formatting
    with open("functional_groups_out.json", "w") as f:
        json.dump(sorted_groups, f, indent=4, sort_keys=True)
    
    print(f"Generated functional_groups_out.json with {len(sorted_groups)} functional group definitions")

# Generate the functional groups output file
generate_functional_groups_output(reactions)


if len(sys.argv) > 1:
    # If a specific reaction is provided as an argument, filter the reactions
    specific_reaction = sys.argv[1]
    
    # Check if the argument is a range (e.g., "60-100")
    if '-' in specific_reaction and not specific_reaction.startswith('-'):
        try:
            start_str, end_str = specific_reaction.split('-', 1)
            start_num = int(start_str)
            end_num = int(end_str)
            
            if start_num > end_num:
                print(f"Invalid range: start ({start_num}) is greater than end ({end_num})")
                sys.exit(1)
            
            # Find reactions within the RXN_NUM range (inclusive)
            reactions = {k: v for k, v in reactions.items() 
                        if v.get("RXN_NUM") is not None and start_num <= v.get("RXN_NUM") <= end_num}
            
            if not reactions:
                print(f"No reactions found with RXN_NUM in range: {start_num}-{end_num}")
                sys.exit(1)
            else:
                found_nums = sorted([v.get("RXN_NUM") for v in reactions.values()])
                print(f"Found {len(reactions)} reactions with RXN_NUM in range {start_num}-{end_num}: {found_nums}")
                
        except ValueError:
            print(f"Invalid range format: '{specific_reaction}'. Use format like '60-100'")
            sys.exit(1)
    
    # Check if the argument is a single number
    elif specific_reaction.isdigit():
        # Find reaction by RXN_NUM
        target_rxn_num = int(specific_reaction)
        reactions = {k: v for k, v in reactions.items() if v.get("RXN_NUM") == target_rxn_num}
        if not reactions:
            print(f"No reaction found with RXN_NUM: {target_rxn_num}")
            sys.exit(1)
    else:
        # Find reaction by name (current behavior)
        reactions = {k: v for k, v in reactions.items() if k == specific_reaction}
        if not reactions:
            print(f"No reaction found with name: {specific_reaction}")
            sys.exit(1)

for reaction_idx, reaction_name in enumerate(reactions):
    if "SKIP" in reaction_name:
        print(f"Skipping reaction: {reaction_name} >> ")
        continue

    reaction_info = reactions[reaction_name]
    functional_groups = reaction_info.get("functional_groups", [])
    if not functional_groups:
        print(f"No functional groups defined for reaction: {reaction_name}")
        continue

    debug_prereactions = reaction_info.get("debug_prereaction", None)
    if debug_prereactions:
        print(f"Found debug prereactions for {reaction_name}: {debug_prereactions}")

    # --- REVISED CALL ---
    reactant_groups = get_forward_reactant_groups(
        functional_groups, reaction_name, debug_prereactions
    )
    num_passed = 0
    num_tested = 0
    num_forward_failures = 0
    num_invalid_reactants = 0  # Track invalid reactants
    num_reverse_skipped = 0   # NEW: Track reverse reactions skipped due to exceptions
    
    # Check if reactant_groups is empty
    if not reactant_groups or not reactant_groups[0]:
        print(f"Warning: No valid reactant molecules found for functional groups in {reaction_name}. Skipping tests.")
        continue

    for reactant_idx, forward_reactants_raw in enumerate(reactant_groups):
        # print(f"\n===== Testing Set {reactant_idx} for {reaction_name} =====")
        
        # Clean reactants first and create Mol objects
        cleaned_reactants = [clean_up_smiles(r) for r in forward_reactants_raw]
        if None in cleaned_reactants:
            print(f"SKIPPING: Invalid reactants found after cleaning for {reaction_name}: {forward_reactants_raw}")
            continue
        
        mols = [Chem.MolFromSmiles(s) for s in cleaned_reactants]
        if any(m is None for m in mols):
            print(f"SKIPPING: Could not create molecules for cleaned reactants: {cleaned_reactants}")
            continue

        # Sanitize and generate both aromatic and Kekule forms
        for mol in mols:
            Chem.SanitizeMol(mol)
        
        aromatic_reactants = [Chem.MolToSmiles(m, canonical=True, kekuleSmiles=False) for m in mols]
        kekule_reactants = [Chem.MolToSmiles(m, canonical=True, kekuleSmiles=True) for m in mols]

        # --- Test Aromatic Form ---
        aromatic_status = execute_and_verify_reaction(
            reaction_name, reaction_info, aromatic_reactants, reactant_idx, "aromatic"
        )
        
        # --- Test Kekule Form ---
        kekule_status = ReactionStatus.SUCCESS
        if aromatic_reactants != kekule_reactants:
            kekule_status = execute_and_verify_reaction(
                reaction_name, reaction_info, kekule_reactants, reactant_idx, "kekule"
            )
        else:
            #  print(f"Kekule form is identical to aromatic for set {reactant_idx}, skipping separate test.")
             kekule_status = aromatic_status
        
        # NEW: Handle different status outcomes including REVERSE_SKIPPED
        if aromatic_status == ReactionStatus.INVALID_REACTANTS or kekule_status == ReactionStatus.INVALID_REACTANTS:
            # Don't count invalid reactants in our test statistics
            num_invalid_reactants += 1
            continue
        
        # Count as tested if we have valid reactants (regardless of reverse skip)
        num_tested += 1
        
        if aromatic_status == ReactionStatus.REVERSE_FAILURE or kekule_status == ReactionStatus.REVERSE_FAILURE:
            # The detailed debug prompt has already been printed.
            raise AssertionError(f"FATAL: Reverse reaction for '{reaction_name}' failed. See debug prompt above.")
        
        # Tally results for the set
        if ((aromatic_status == ReactionStatus.SUCCESS or aromatic_status == ReactionStatus.REVERSE_SKIPPED) and 
            (kekule_status == ReactionStatus.SUCCESS or kekule_status == ReactionStatus.REVERSE_SKIPPED)):
            num_passed += 1
            # Count reverse skipped reactions for reporting
            if aromatic_status == ReactionStatus.REVERSE_SKIPPED or kekule_status == ReactionStatus.REVERSE_SKIPPED:
                num_reverse_skipped += 1
        elif aromatic_status == ReactionStatus.FORWARD_FAILURE or kekule_status == ReactionStatus.FORWARD_FAILURE:
            num_forward_failures += 1
            print(f"NOTE: Set {reactant_idx} for {reaction_name} logged as a FORWARD FAILURE.")

    # Final summary for the reaction
    print(f"\nFINAL SUMMARY for {reaction_name}:")
    print(f"  - Total Sets Tested: {num_tested}")
    print(f"  - Passed Sets: {num_passed}")
    print(f"  - Forward Reaction Failures: {num_forward_failures}")
    print(f"  - Invalid Reactants (Skipped): {num_invalid_reactants}")
    print(f"  - Reverse Reactions Skipped (Exception): {num_reverse_skipped}")  # NEW: Report skipped reverse reactions

    if num_passed == 0 and num_tested > 0:
        print(f"\nWARNING: No reactant sets fully passed validation for {reaction_name}.")
    elif num_tested > 0:
        print(f"\nSUCCESS: {num_passed}/{num_tested} reactant sets passed for {reaction_name}.\n")
    
    # NEW: Additional reporting for reverse skipped reactions
    if num_reverse_skipped > 0:
        print(f"NOTE: {num_reverse_skipped} reactant sets had reverse reactions skipped due to exception pattern.\n")