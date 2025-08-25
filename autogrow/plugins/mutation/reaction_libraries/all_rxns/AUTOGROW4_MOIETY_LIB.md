# Fragment Library Generation Notes

## Source and Initial Filtering

The fragment libraries included with AutoGrow4 were derived from a subset of the
ZINC15 database, retrieved on December 19, 2019.

The initial query from ZINC15 was constrained to compounds meeting the following
criteria:

* **Commercial Availability**: Readily available for purchase ("Wait OK" status
  or sooner).
* **Molecular Weight (MW)**: Less than or equal to 250 Da.
* **LogP**: Less than or equal to 5.0.
* **Reactivity**: All reactivity levels were included (e.g., "Hot").

This query returned an initial set of 19,274,338 substances.

## Secondary Filtering and Curation

This initial list of compounds was subjected to further computational filtering:

1. Each SMILES string was processed to ensure it could be successfully parsed
   into a valid RDKit molecule object.
2. The molecules were then passed through a strict Lipinski filter to select for
   drug-like properties.

The resulting high-quality molecules were used to populate the initial
functional group libraries. The compounds within each library were sorted by
molecular weight. To maintain a focus on small, extensible fragments, any
library containing more than 5,000 molecules was truncated to include only the
5,000 compounds with the lowest molecular weight. This preference for smaller
starting fragments is intentional, as the AutoGrow process is designed to
increase molecular size.

Finally, all compounds were validated to confirm their ability to participate in
the chemical reactions for which their respective moiety libraries are used.

## Final Library Composition

The `all_rxns` library was created by merging the complementary molecule
libraries from the `robust_rxns` and `click_chem_rxns` sets. After removing
redundancies, the final combined set consists of **185,615** unique compounds.

*The default fragment libraries included with AutoGrow were derived from a
subset of the ZINC database (https://zinc.docking.org/). We thank ZINC for
allowing us to distribute these fragment libraries to AutoGrow users.*