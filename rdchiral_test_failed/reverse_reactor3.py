import json
from typing import List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
import re
import gzip
import random
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import os
import sys

# Assume reactions.json exists and is in the correct format
with open("reactions.json", "r") as file:
    contents = file.read()
    # Replace some placeholders
    contents = contents.replace("AZIDE_SEARCH", "[$(N=[N+]=[N-]),$([N-][N+]#N)]")
    reactions = json.loads(contents)

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

    # Strange tetrazole resonance structure I see in some SMILES
    smiles = smiles.replace("C2=N[N+]=NN2", "C2=NN=NN2")
    # Add more replacements if needed

    # Also, smiles should be canonicalized
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid SMILES: ", smiles)
        return None
    
    try:
        # Let RDKit handle hydrogen assignment properly
        Chem.SanitizeMol(mol)
        # Remove explicit hydrogens and let RDKit recalculate
        mol = Chem.RemoveHs(mol)
        smiles = Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        print(f"Sanitization failed for {smiles}: {e}")
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

# --- MODIFIED FUNCTION ---
def run_rxn(rxn_string: str, reactants_smi_list: List[str])-> List[List[str]]:
    """
    Runs a reaction with a list of reactant SMILES and returns the product sets.
    
    Args:
        rxn_string: The SMARTS reaction string.
        reactants_smi_list: A list of SMILES strings for the reactants.
        
    Returns:
        A list of product sets, where each product set is a list of SMILES strings.
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
                # Sanitize the molecule to avoid errors in MolToSmiles
                # Chem.SanitizeMol(product)
                smiles = Chem.MolToSmiles(product)
                smiles = clean_up_smiles(smiles)
                if smiles:
                    product_smiles.append(smiles)
                else:
                    print(f"Warning: Failed to clean up product SMILES: {Chem.MolToSmiles(product)}")
                    print(f"Reaction string: {rxn_string}")
                    print(f"Reactants: {reactants_smi_list}")
                    print("")
        if product_smiles:
            # Sort the smiles in each set to make comparisons order-independent
            predicted_products.append(sorted(product_smiles))
            
    return predicted_products

def apply_prereaction(smiles: str, prereaction_strings: List[str]) -> str:
    """
    Apply a list of prereaction strings to standardize a molecule.
    If any prereaction yields products, return the first product from the first successful reaction.
    If no prereaction yields products, return the original SMILES.
    
    Args:
        smiles: The input SMILES string
        prereaction_strings: List of SMARTS reaction strings to try
        
    Returns:
        Standardized SMILES string
    """
    for prereaction_string in prereaction_strings:
        try:
            # Try to run this prereaction on the molecule
            products = run_rxn(prereaction_string, [smiles])
            
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
    print(f"  No prereaction applied to: {smiles}")
    return smiles

def get_forward_reactant_groups(functional_groups: List[str], debug_prereactions: Optional[List[str]] = None) -> List[List[str]]:
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
    for i in range(len(reactant_groups[0])):
        # First, do the forward reaction
        forward_reactants = [clean_up_smiles(reactant_groups[j][i]) for j in range(len(reactant_groups))]
        forward_reactants = [f for f in forward_reactants if f]
        
        # Apply debug prereactions if they exist
        if debug_prereactions:
            print(f"Applying debug prereactions to reactant set {i}:")
            standardized_reactants = []
            for j, reactant in enumerate(forward_reactants):
                # Apply prereactions to each reactant individually
                standardized_reactant = apply_prereaction(reactant, debug_prereactions)
                standardized_reactants.append(standardized_reactant)
            forward_reactants = standardized_reactants
        
        reactant_groups2.append(forward_reactants)

    return reactant_groups2

def save_reaction_visualization(reaction_name: str, forward_reactants: List[str], forward_products: List[List[str]], 
                               all_reverse_products: List[str], original_reactants_set: set, 
                               reverse_products_set: set, reactant_index: int, 
                               success: bool = True) -> None:
    """
    Save a visualization of the reaction pathway as PNG.
    
    Args:
        reaction_name: Name of the reaction
        forward_reactants: List of reactant SMILES
        forward_products: List of product sets from forward reaction
        all_reverse_products: List of reverse reaction products
        original_reactants_set: Set of cleaned original reactant SMILES
        reverse_products_set: Set of cleaned reverse product SMILES
        reactant_index: Index of the current reactant set
        success: Whether the reaction was successful
    """
    # Create a directory for reaction visualizations if it doesn't exist
    viz_dir = "reaction_visualizations"
    os.makedirs(viz_dir, exist_ok=True)

    # Create a safe filename
    safe_reaction_name = re.sub(r'[^\w\-_\.]', '_', reaction_name)
    
    # Get actual molecules for original reactants
    reactant_mols = [Chem.MolFromSmiles(smi) for smi in forward_reactants if Chem.MolFromSmiles(smi) is not None]
    
    # Get product molecules from forward reaction
    product_mols = []
    for reverse_reactants_set in forward_products:
        for product_smi in reverse_reactants_set:
            mol = Chem.MolFromSmiles(product_smi)
            if mol is not None:
                product_mols.append(mol)

    # Get reverse product molecules - use the cleaned set that's actually being compared
    reverse_mols = []
    for smi in reverse_products_set:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            reverse_mols.append(mol)

    # Create comprehensive visualization showing the full reaction pathway
    all_mols = reactant_mols + [None] + product_mols + [None] + reverse_mols
    legends = ([f"Reactant {i+1}" for i in range(len(reactant_mols))] + 
              ["→"] + 
              [f"Product {i+1}" for i in range(len(product_mols))] + 
              ["Reverse→"] + 
              [f"Reverse {i+1}" for i in range(len(reverse_mols))])

    # Add information about missing reactants if failed
    if not success:
        missing_reactants = original_reactants_set - reverse_products_set
        print(f"Missing from reverse products: {missing_reactants}")

    # Determine the status for filename
    status = "PASS" if success else "FAIL"
    
    # Save the comprehensive reaction pathway as PNG
    png_filename = f"{viz_dir}/{safe_reaction_name}_{reactant_index:03d}_{status}.png"
    
    try:
        # RDKit's MolsToGridImage returns a PIL Image when useSVG=False
        img = Draw.MolsToGridImage(all_mols, molsPerRow=max(len(reactant_mols), len(product_mols), len(reverse_mols)) + 1, 
                                  subImgSize=(200, 200), legends=legends, useSVG=False)
        # Save the PIL Image directly
        img.save(png_filename, "PNG")
        # print(f"Saved reaction visualization: {png_filename}")
    except Exception as e:
        print(f"Failed to save PNG visualization {png_filename}: {e}")
        # Fallback to SVG if PNG fails
        try:
            svg_filename = f"{viz_dir}/{safe_reaction_name}_{reactant_index:03d}_{status}.svg"
            img_svg = Draw.MolsToGridImage(all_mols, molsPerRow=max(len(reactant_mols), len(product_mols), len(reverse_mols)) + 1, 
                                          subImgSize=(200, 200), legends=legends, useSVG=True)
            with open(svg_filename, 'w') as f:
                f.write(img_svg)
            print(f"Saved SVG fallback: {svg_filename}")
        except Exception as e2:
            print(f"Failed to save SVG fallback: {e2}")

if len(sys.argv) > 1:
    # If a specific reaction is provided as an argument, filter the reactions
    specific_reaction = sys.argv[1]
    
    # Check if the argument is a number
    if specific_reaction.isdigit():
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
    # if reaction_idx < 84:
    #     continue
    if "SKIP" in reaction_name:
        print(f"Skipping reaction: {reaction_name} >> ")
        continue

    # Get the reaction information
    reaction_info = reactions[reaction_name]

    functional_groups = reaction_info.get("functional_groups", [])
    if not functional_groups:
        print(f"No functional groups defined for reaction: {reaction_name}")
        continue

    # Get debug prereactions if they exist
    debug_prereactions = reaction_info.get("debug_prereaction", None)
    if debug_prereactions:
        print(f"Found debug prereactions for {reaction_name}: {debug_prereactions}")

    reactant_groups = get_forward_reactant_groups(functional_groups, debug_prereactions)

    reaction_string = reaction_info["reaction_string"]
    reverse_reaction_strings = reaction_info.get("reverse_reaction_strings", [])
    num_passed = 0
    
    for reactant_idx, forward_reactants in enumerate(reactant_groups):
        # First, do the forward reaction
        forward_reactants = [clean_up_smiles(r) for r in forward_reactants]
        if None in forward_reactants:
            print(f"Invalid reactants found in forward reaction for {reaction_name}: {forward_reactants}")
            continue
        forward_products = run_rxn(reaction_string, forward_reactants)

        # Make sure none of the products are None
        forward_product_is_none = False
        for t1 in forward_products:
            for t2 in t1:
                if t2 is None:
                    forward_product_is_none = True
                    break
            if forward_product_is_none:
                break
        if forward_product_is_none:
            print(f"Invalid products found in forward reaction for {reaction_name}: {forward_products}")
            continue

        if len(forward_products) == 0:
            print(f"No products generated for forward reaction: {reaction_string} with reactants: {forward_reactants}")
            # assert False, f"No products generated for forward reaction: {reaction_string} with reactants: {forward_reactants}"
            continue

        # Now, apply each of the reverse reactions to each of the products, collecting them all in a single set
        all_reverse_products = []
        for reverse_reactants_set in forward_products:
            for reverse_reaction_string in reverse_reaction_strings:
                # For each product in the current product set, run the reverse reaction
                for reverse_reactant in reverse_reactants_set:
                    reverse_products = run_rxn(reverse_reaction_string, [reverse_reactant])
                    for r in reverse_products:
                        all_reverse_products.extend(iter(r))
        # Convert everything to sets
        all_reverse_products_before_cleaning = set(all_reverse_products[:])
        original_reactants_set = {clean_up_smiles(p) for p in forward_reactants if clean_up_smiles(p)}
        reverse_products_set = {clean_up_smiles(p) for p in all_reverse_products if clean_up_smiles(p)}

        if original_reactants_set.issubset(reverse_products_set):
            print("Passed!") # , reaction_name, reverse_products_set)
            num_passed += 1
            
            # Save visualization for successful reaction
            save_reaction_visualization(reaction_name, forward_reactants, forward_products, 
                                      all_reverse_products, original_reactants_set, 
                                      reverse_products_set, reactant_idx, success=True)
        else:
            # Save visualization for failed reaction
            save_reaction_visualization(reaction_name, forward_reactants, forward_products, 
                                      all_reverse_products, original_reactants_set, 
                                      reverse_products_set, reactant_idx, success=False)
            
            print(f"Missing reactants: {original_reactants_set - reverse_products_set}")

            # Some of the original reactants were not found in the reverse products. The reverse reaction failed.
            print(f"Reverse reaction for '{reaction_name}' FAILED.\n"
                f"  > Original Reactants: {original_reactants_set}\n"
                f"  > Reverse Products: {reverse_products_set}")
            print(f"  > Reverse Reaction Strings: {reverse_reaction_strings}")
            print(f"  > Forward Reaction String: {reaction_string}")
            print(f"  > Reactants: {forward_reactants}")
            print(f"  > Products: {forward_products}")
            print("")
            print("")
            print("Consider this forward reaction:")
            print(reaction_string)
            print("With these reactants:")
            print(forward_reactants)
            print("I get this products:")
            print({item for sublist in forward_products for item in sublist})
            print("I then take that product and feed it through these reverse reactions:")
            print(reverse_reaction_strings)
            print("And I get these products:")
            print(set(all_reverse_products_before_cleaning))
            print("I would hope that the original reactants are a subset of the reverse products, but they are not.")
            print("These original reactants were not found in the reverse products:")
            print(original_reactants_set - reverse_products_set)
            print("This is a failure of the reverse reaction. Can you fix it?")
            print("")
            raise AssertionError(f"Reverse reaction for '{reaction_name}' FAILED. Check the output above for details.")

    if num_passed == 0:
        assert False, f"Warning: No successful reactions for {reaction_name}. This might indicate an issue with the reaction setup or the reactants."