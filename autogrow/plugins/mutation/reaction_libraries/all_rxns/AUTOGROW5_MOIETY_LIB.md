# Moiety Database Generation Pipeline

This README.md file describes how the AutoGrow5 moiety database was generated.
The database contains multiple gzipped SMILES files, where each file corresponds
to a specific chemical moiety (i.e., a functional group). These files contain
curated lists of unique, small-molecule compounds that contain the corresponding
moiety.

The generation process involved a multi-step pipeline that collected, filtered,
categorized, and augmented molecules to create the final libraries.

## Source Compounds

The first step involves gathering a large, diverse set
of initial candidate molecules from two sources.

* **NCI/DTP Open Chemicals Repository**: The process begins by automatically
  downloading the "Chem2D_Jun2016" compound collection from the NCI/DTP Open
  Chemicals Repository. This dataset provides a substantial base of over 280,000
  compounds. Per the NCI website, these compounds are [in the public
  domain](https://dctd.cancer.gov/drug-discovery-development/reagents-materials/vialed-plated-compounds-0).
* **AutoGrow4 Libraries**: To enrich the initial set, the pipeline also fetches
  the original AutoGrow4 moiety database from the AutoGrow4 repository. This
  database is derived from a small subset of the [ZINC
  database](https://zinc.docking.org/), used by permission granted on December
  13, 2019. We thank ZINC for allowing us to distribute these fragment libraries
  to AutoGrow users.

During this initial collection, a first-pass filter is applied to all incoming molecules:

1. Molecules are "desalted" by isolating their largest covalent fragment.
2. Fragments with a molecular weight greater than 250.0 Da are discarded.
3. A canonical SMILES string is generated for each remaining fragment to ensure
   only unique molecules are retained.

## Compound Filtering

After the initial collection, a second, more stringent set of filters is
applied. Molecules that match any of the following criteria are discarded:

* Four or more rotatable bonds
* cLogP greater than 3.0
* Compounds with atoms that are Ac, Ag, Al, Am, As, Au, Ba, Be, Bi, Bk, Ca, Cd,
  Ce, Cf, Cm, Co, Cr, Cs, Cu, Dy, Er, Es, Eu, Fe, Fm, Fr, Ga, Gd, Ge, Hf, Hg,
  Ho, In, Ir, K, La, Li, Lr, Lu, Md, Mg, Mn, Mo, Na, Nb, Nd, Ni, No, Np, Os, Pa,
  Pb, Pd, Pm, Po, Pr, Pt, Pu, Ra, Rb, Re, Rh, Ru, Sb, Sc, Se, Si, Sm, Sn, Sr,
  Ta, Tb, Tc, Te, Th, Ti, Tl, Tm, U, V, W, Y, Yb, Zn, Zr
* Compounds that RDKit cannot parse or sanitize
* Compounds that RDKit cannot neutralize (i.e., adding or removing hydrogens to
  neutralize formal charges on atoms)

## Moiety Categorization

Next, each fully-filtered compound is categorized based on the functional groups
it contains. This categorization is formed by matching each molecule against the
moiety definitions (SMARTS strings) specified in the `reactions.json` file.

Since a single molecule can contain multiple functional groups, it may be
assigned to several moiety categories.

## Library Augmentation via Chemical Synthesis and Retrosynthesis

The final step addresses moiety libraries containing fewer than 5,000 molecules.
These underpopulated libraries are expanded using an iterative, *in silico*
synthesis engine that attempts to generate new, valid molecules using two
strategies:

* **Elaboration (Forward Synthesis)**: A small, existing molecule from an
  underpopulated library is chosen as a starting point. A second "partner"
  molecule from another library is selected to react with it, guided by the
  reaction definitions in `reactions.json`. If the resulting product still
  contains the target moiety and passes all filters above, it is added to the
  database.
* **Decomposition (Retrosynthesis)**: A larger molecule from a library is
  selected, and the retrosynthetic rules from `reactions.json` are applied to
  break it down into smaller precursor molecules. If any of these smaller
  precursors are valid new molecules that still contain the target moiety, they
  are added to the appropriate libraries.

This augmentation loop continues until all moiety libraries meet the target
count or until the system can no longer generate new valid molecules. In
practice, it is not possible to generate more than 5,000 molecules for every
moiety while adhering to the strict filtering criteria. Each moiety library is
saved to a gzipped SMILES file.
