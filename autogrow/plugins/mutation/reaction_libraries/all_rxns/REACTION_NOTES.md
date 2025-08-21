# Reaction Library Notes

This document details the origin and modifications of the AutoGrow4 reaction
sets, which were subsequently adapted for use in AutoGrow5. The final library is
a merger of two distinct collections: the "Click Chemistry" set and the "Robust"
reaction set.

## Click Chemistry Reaction Set (`click_chem_rxns`)

This reaction set is based on the reactions used in AutoClickChem, which were
inspired by Click Chemistry principles. All reaction definitions were converted
from their original 3D format into 2D SMARTS representations for broader
compatibility.

Several key modifications were made to improve chemical accuracy and generality:

1. **Corrected Moiety Definitions for Carboxyl and Thiol Groups**:

* In a previous implementation, carboxyl groups (`-COOH`) were sometimes
  misidentified as alcohols, leading to incorrect reactions. This has been
  corrected; carboxyl groups are no longer treated as simple alcohols.
* Similarly, thiol-carboxyl groups (`-COSH`) are no longer misidentified as
  standard thiols or alcohols.

2. **Consolidated Halide Reactions**:

* AutoClickChem previously used separate reactions for acyl, alkyl, and aryl
  halides because their underlying chemical mechanisms (e.g., Nucleophilic Acyl
  Substitution vs. SN1/SN2 vs. SNAr) are different.
* However, from a structural transformation perspective represented by SMARTS,
  these reactions are very similar. For example, the conversion of a halide to
  an azide (`R-X` → `R-N₃`) can be described with a single, more general SMARTS
  pattern that covers multiple halide types.
* We have therefore condensed reactions for different halide types into a
  single, more general reaction definition. This applies to the following
  transformations:
  * Halide → Cyanide
  * Halide → Azide
  * Halide + Thiol → Thioether
  * Halide + Alcohol → Ether
  * Halide + Amine → Amine

3. **Refined Reactivity for SP2-Hybridized Halides**:

* The previous version treated any halide attached to an sp2-hybridized carbon
  (including vinyl halides) as an acyl halide.
* The modern implementation correctly excludes vinyl halides, which are poor
  substrates for substitution reactions, while still allowing reactions with
  aromatic and acyl halides.

### Additional Notes

**Non-Cyclic Acid Anhydrides**: The reactions involving acid anhydrides are
restricted to non-cyclic variants. While cyclic anhydrides can react in a wet
lab setting, the current SMARTS definitions are not equipped to handle the
ring-breaking transformations that would occur, as this would alter the
molecular graph in a way that is difficult to represent. This limitation is
consistent with previous versions of AutoGrow and AutoClickChem.

## Robust Reaction Set (`robust_rxns`)

This set of reactions is derived from the work of Hartenfeller et al. and has
undergone several modifications to improve its utility and robustness.

*Hartenfeller, M. et al. J. Chem. Inf. Model. 2011, 51 (12), 3093-3098.*

*Hartenfeller, M. et al. J. Chem. Inf. Model. 2012, 52 (5), 1167-1178.*

### Minor Changes

* SMARTS patterns were adjusted to be more inclusive where chemically
  appropriate.
* Functional group definitions were canonicalized. For instance, multiple
  reactions that required a general alkyne now use a single, standardized alkyne
  SMARTS pattern, allowing them to draw from the same complementary moiety
  library.
* Minor errors in specific reaction examples and definitions were corrected
  (e.g., typos in Reactions #22 and #53).

### Major Changes

1. **Reaction #48 (Sulfonamide Synthesis)**:

* The original reaction used sulfonyl chloride as a reactant. However, sulfonyl
  chloride is often an unstable intermediate that is not commercially available.
* In a typical laboratory synthesis, the reaction starts with a more stable and
  available sulfonic acid, which is converted to sulfonyl chloride *in situ*
  before reacting with an amine.
* To better reflect the availability of starting materials for computational
  generation, we have simplified this two-step process into a single
  computational reaction. The reaction now starts directly with a **sulfonic
  acid** and an amine to produce the sulfonamide.

2. **Reaction #7 (Thiazole Synthesis)**:

* The original SMARTS definition for the required thioamide reactant was highly
  restrictive and matched compounds that are relatively rare in commercial
  databases.
* Thioamides are easily synthesized from more common amides. Therefore, the
  reaction was expanded to accept both **amides and thioamides** as valid
  starting materials.
* Furthermore, the SMARTS definition was updated to recognize both tautomeric
  forms of the thioamide reactant (`NH₂-C(R)=S` and `NH=C(R)-SH`), making the
  pattern more robust.

## Combined Reaction Library

The final directory of reactions is a merger of the **Click Chemistry** and
**Robust** libraries.

To help identify the origin of each functional group, the following naming
convention is used for their corresponding moiety libraries:

* Groups from the Robust set end with `_robust`.
* Groups from the Click Chemistry set end with `_clickchem`.