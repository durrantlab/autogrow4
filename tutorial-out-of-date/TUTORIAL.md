# AutoGrow5 Tutorial: Evolutionary Drug Design

Welcome to this comprehensive tutorial on AutoGrow5, a powerful tool for computational drug discovery that uses evolutionary algorithms to generate novel molecules with optimized properties. In this tutorial, you'll learn how to install, configure, and run AutoGrow5 to design molecules tailored to your specific protein targets.

## What is AutoGrow5?

AutoGrow5 is an open-source program that mimics biological evolution to discover new drug-like molecules. Starting with a collection of seed compounds, it iteratively generates new molecular variants through chemical modifications (mutations) and combinations (crossovers), then evaluates each generation's fitness through molecular docking. Over multiple generations, this process evolves increasingly optimized compounds for your target protein.

## Before We Begin: Understanding File Paths

Throughout this tutorial, you'll encounter placeholder paths marked with `ALL_CAPS` formatting. These must be replaced with actual paths on your system. For example:

```bash
cd /PATH/TO/DESIRED_DIR/
```

Should become something like:

```bash
cd /home/username/Documents/AutoGrow5/
```

For simplicity, we'll refer to the main AutoGrow5 directory as `/autogrow5/` throughout this tutorial, but remember to use your actual installation path.

## System Requirements

Before installing AutoGrow5, ensure your system meets these requirements:

- **Ubuntu 16.04 or higher** (recommended and fully supported)
- **macOS 10.13 High Sierra or higher** (tested and supported)
- **Windows** (may work but is not officially supported)

## Step 1: Installing AutoGrow5

Let's start by setting up AutoGrow5 and its dependencies. The easiest approach uses conda to create an isolated environment with all necessary components.

### Setting Up the Conda Environment

AutoGrow5 comes with an `environment.yml` file that automatically installs Python 3.7+ and all required Python libraries. This includes essential packages like:

- **RDKit**: The core cheminformatics library for molecular manipulation
- **NumPy, pandas, scikit-learn**: For data processing and machine learning
- **Matplotlib**: For visualizing results
- **ProLIF**: For analyzing protein-ligand interactions

To create and activate the environment:

```bash
conda env create -f environment.yml
conda activate AutoGrow5
```

This single command handles most of the installation process, but we still need to install a few external dependencies.

## Step 2: Installing Required External Dependencies

While conda handles the Python packages, AutoGrow5 relies on two external programs that must be installed separately.

### Installing Open Babel

Open Babel is essential for converting between different molecular file formats, particularly for transforming 2D SMILES strings into 3D molecular structures needed for docking.

**Installation varies by operating system:**

- **Ubuntu/Debian users**: `sudo apt-get install openbabel`
- **macOS users with Homebrew**: `brew install openbabel`
- **Other systems**: Follow the installation guide at https://openbabel.org/docs/dev/Installation/install.html

**Finding Open Babel's location:**

After installation, you'll need to tell AutoGrow5 where to find Open Babel:

```bash
which obabel
```

Remember this path—you'll use it later with the `--obabel_path` parameter.

### Installing Docking Software

AutoGrow5 evaluates molecular fitness through docking simulations, so you need to install at least one supported docking program. You have several excellent options:

**AutoDock Vina Family (Recommended for beginners):**

- **AutoDock Vina**: The original, widely-used program. Download from [The Scripps Research Institute](http://vina.scripps.edu/download.html)
- **QuickVina 2**: A faster variant ideal for large-scale screening. Get it from [GitHub](https://github.com/QVina/qvina/releases)
- **Smina**: Another optimized fork with additional features. Available on [SourceForge](https://sourceforge.net/projects/smina/files/)

**GNINA (Advanced option):**

- **GNINA**: A deep learning-enhanced docking program that often provides superior results. Find installation instructions on the [GNINA GitHub repository](https://github.com/gnina/gnina)

Choose the docking program that best fits your needs and computational resources. Note the installation path—you'll specify this with the `--docking_executable` parameter.

## Step 3: Understanding How AutoGrow5 Works

Before running your first experiment, let's explore how AutoGrow5's evolutionary algorithm operates. Understanding this process will help you configure parameters effectively.

### The Evolutionary Cycle

AutoGrow5 evolves molecules through repeated cycles of selection, breeding, and evaluation:

**1. Selection Phase**
The algorithm selects promising molecules from the current generation to serve as "parents." This selection uses `Selector` plugins that can prioritize molecules based on docking scores, structural diversity, or custom criteria.

**2. Breeding Phase**
Selected parents generate new "offspring" molecules through two mechanisms:

- **Mutation**: `Mutation` plugins create variants by making small chemical modifications to single parent molecules
- **Crossover**: `Crossover` plugins combine fragments from two different parent molecules

**3. Elitism**
The best-performing molecules from the previous generation automatically advance to ensure promising candidates aren't lost.

**4. Quality Control**
All new molecules pass through `SmilesFilter` plugins that remove compounds with undesirable properties like poor drug-likeness or reactive functional groups.

**5. 3D Structure Generation**
Molecules that pass filtering are converted from 2D SMILES notation to 3D structures using `SmiTo3DSdf` plugins (typically Open Babel).

**6. Docking and Evaluation**
The 3D structures undergo molecular docking with your target protein using `Docking` plugins. Additional `PoseFilter` and `Rescoring` plugins can refine the evaluation process.

**7. Ranking and Selection**
All molecules in the new generation are ranked by fitness score, creating the foundation for the next evolutionary cycle.

This process repeats for your specified number of generations, progressively optimizing the molecular population.

### Plugin Architecture

AutoGrow5's modular design makes it highly customizable. Each step in the workflow is handled by specific plugin types:

**Docking Plugins:**

- `VinaLikeDocking`: Interfaces with AutoDock Vina and compatible programs
- `GNINADocking`: Integrates GNINA for deep learning-enhanced docking

**Scoring Enhancement:**

- `LigandEfficiency`: Calculates efficiency-based scores (docking score ÷ heavy atom count)

**Molecular Conversion:**

- `ObabelSmiTo3DSDF`: Converts SMILES to 3D SDF format using Open Babel

**Quality Filters:**

- `LipinskiStrictFilter`: Enforces Lipinski's Rule of Five for drug-likeness
- `PAINSFilter`: Removes Pan-Assay Interference Compounds
- `MolWtFilter`: Filters by molecular weight ranges

**Selection Strategies:**

- `RankSelector`: Chooses top-scoring molecules
- `RouletteSelector`: Uses probabilistic selection based on fitness

**Molecular Evolution:**

- `FragmentAddition`: Adds chemical fragments using reaction libraries
- `FragmentRemoval`: Removes fragments through reverse reactions
- `MergeMCS`: Combines molecules based on maximum common substructures

**Interaction Analysis:**

- `PoseFilter`: Evaluates specific protein-ligand interactions

**Parallel Processing:**

- `ShellParallelizer`: Manages concurrent execution of external programs

## Step 4: Preparing Your Protein Target

Successful molecular evolution requires a properly prepared protein structure. Follow these critical steps:

### Cleaning Your PDB File

Start with a high-quality crystal structure and remove unnecessary components:

1. **Remove non-protein atoms**: Delete ligands, water molecules, and ions using molecular visualization software like PyMOL or VMD
2. **Select relevant chains**: For multi-chain proteins, keep only chains containing or near your target binding site
3. **Add hydrogen atoms**: Crystal structures lack hydrogens. Use tools like [PDB2PQR](http://nbcr-222.ucsd.edu/pdb2pqr_2.0.0/) to add them at physiological pH
4. **Convert file format**: If PDB2PQR outputs a PQR file, convert it back to PDB format using Open Babel

### Defining the Binding Site

Docking requires a defined search space around your target binding site. You need six parameters:

- `--center_x`, `--center_y`, `--center_z`: The geometric center of your binding site
- `--size_x`, `--size_y`, `--size_z`: The dimensions of the search box

**Tip**: The Python library `scoria` can help calculate these coordinates from selected binding site residues.

## Step 5: Running Your First AutoGrow5 Experiment

Now let's run AutoGrow5 with a practical example. You can provide parameters either through command-line arguments or a JSON configuration file.

### Command-Line Execution

Here's a complete example command:

```bash
python run_autogrow.py \
    --receptor_path /PATH/TO/receptor.pdb \
    --center_x -70.76 --center_y 21.82 --center_z 28.33 \
    --size_x 25.0 --size_y 16.0 --size_z 25.0 \
    --source_compound_file /PATH/TO/source_compounds.smi \
    --output_directory /PATH/TO/output_directory/ \
    --num_generations 10 \
    --number_of_mutants 50 \
    --number_of_crossovers 50 \
    --top_mols_to_seed_next_generation 25 \
    --number_elitism_advance_from_previous_gen 10 \
    --procs_per_node -1 \
    --VinaLikeDocking \
    --docking_executable /PATH/TO/vina \
    --obabel_path /PATH/TO/obabel \
    --LipinskiStrictFilter
```

### Using JSON Configuration Files

For complex experiments, JSON files offer better organization:

```bash
python run_autogrow.py --json /PATH/TO/params.json
```

### Understanding Key Parameters

Let's break down the most important parameters that control population dynamics:

**Population Size Control:**

- `--number_of_mutants`: New molecules created via mutation each generation
- `--number_of_crossovers`: New molecules created via crossover each generation  
- `--number_elitism_advance_from_previous_gen`: Elite molecules carried over unchanged

The total population equals the sum of these three values.

**Parent Selection:**

- `--top_mols_to_seed_next_generation`: Number of top-scoring molecules selected as parents
- `--diversity_mols_to_seed_first_generation`: Molecules selected for structural diversity in generation 1
- `--diversity_seed_depreciation_per_gen`: How diversity selection decreases in later generations

**First Generation Adjustments:**

You can use different population sizes for the initial generation with parameters like `--number_of_mutants_first_generation`.

### Getting Help

AutoGrow5 offers extensive customization options. View all available parameters:

```bash
python run_autogrow.py --help
```

## Step 6: Preparing Input Files

### Source Compound Formats

AutoGrow5 accepts starting compounds in two formats:

- **Tab-delimited `.smi` files**: Simple text files with SMILES strings
- **`.sdf` files**: 3D molecular structure files

Example files are available in `/autogrow/source_compounds/` to help you get started.

### Using Pre-evaluated Compounds

To establish a baseline, use the `--process_input_compounds` flag to dock and evaluate your source compounds as "generation 0."

If your `.smi` file already contains docking scores from previous runs, format them as:

```text
SMILES <tab> ID <tab> fitness_score <tab> docking_score ...
```

AutoGrow5 will reuse these scores, saving computation time.

## Step 7: Optimizing Performance with Parallelization

AutoGrow5 can leverage multiple CPU cores to accelerate calculations.

### Python Operations

For internal operations like mutation and crossover:

- `--procs_per_node -1`: Use all available CPU cores
- `--procs_per_node N`: Use N specific cores
- `--multithread_mode serial`: Force single-threaded execution if needed

### External Program Parallelization

For docking and other external programs, choose from these `ShellParallelizer` options:

- **`PythonMultiprocessing`** (default): Standard Python parallelization
- **`ParallelExec`**: Uses GNU Parallel for enhanced performance
- **`Slurm`**: For high-performance computing clusters

**HPC Cluster Example:**

```bash
python run_autogrow.py ... \
    --Slurm \
    --slurm_template_file /PATH/TO/template.sh \
    --sbatch_path /PATH/TO/sbatch
```

## Step 8: Advanced Customization

### Creating Custom Plugins

AutoGrow5's plugin system allows you to extend functionality. To create a custom filter plugin:

1. **Create a Python file** in the appropriate directory (e.g., `autogrow/plugins/smiles_filters/`)
2. **Define a class** inheriting from the base class (e.g., `SmilesFilterBase`)
3. **Implement required methods**:
   - `run_filter()`: Your main filtering logic
   - `add_arguments()`: Command-line flag definition
4. **Enable your plugin** using the flag you defined

AutoGrow5 automatically discovers new plugins, making them immediately available.

### Building Custom Reaction Libraries

For specialized chemistry, you can create custom reaction libraries for the `FragmentAddition` plugin.

**Required Components:**

1. **`rxn_library.json`**: Defines your chemical reactions
2. **`complementary_mols/` directory**: Contains building block molecules

**Creating the JSON file:**

Each reaction needs these specifications:

- `reaction_name`: Unique identifier
- `num_reactants`: Number of reactants (1 or 2)
- `reaction_string`: SMARTS representation of the transformation
- `functional_groups`: Names of required functional groups
- `group_smarts`: SMARTS patterns for each functional group
- `example_rxn_reactants` and `example_rxn_product`: Validation examples
- `RXN_NUM`: Unique integer identifier

**Building Block Files:**

Create compressed SMILES files (`.smi.gz`) for each functional group:

- Filename must match the functional group name
- Format: `SMILES_STRING <tab> MOLECULE_ID <tab> MOLECULAR_WEIGHT`

**Using Your Custom Library:**

```bash
python run_autogrow.py ... \
    --FragmentAddition \
    --rxn_library_path /PATH/TO/YOUR/CUSTOM_LIBRARY/
```

**Testing Your Library:**

Validate your custom reactions before use:

```bash
python autogrow/accessory_scripts/test_complementary_mol_library.py \
    --rxn_library_path /PATH/TO/YOUR/CUSTOM_LIBRARY/ \
    --output_folder /PATH/TO/TEST_OUTPUT/
```

## Conclusion

You now have the knowledge to harness AutoGrow5's evolutionary approach for drug discovery. Start with simple experiments using the default settings, then gradually incorporate advanced features like custom plugins and reaction libraries as your expertise grows.

Remember that successful computational drug discovery combines algorithmic power with chemical intuition. Use AutoGrow5's results as starting points for further optimization, and always validate promising compounds through experimental testing.

Happy molecular evolution!
