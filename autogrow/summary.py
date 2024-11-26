from autogrow.types import Compound
from autogrow.utils.logging import log_info
import numpy as np
import math


def generate_summary_html(output_dir: str) -> None:
    """
    Generate a standalone HTML visualization of AutoGrow4 results.
    
    Args:
        output_dir (str): Path to AutoGrow4 output directory containing generation_X folders
    """
    import os
    import glob
    import json

    # Read all generation data
    all_generations = []
    generation_dirs = sorted(glob.glob(os.path.join(output_dir, "generation_*")))

    for gen_dir in generation_dirs:
        gen_file = os.path.join(gen_dir, f"{os.path.basename(gen_dir)}_ranked.smi")
        if not os.path.exists(gen_file):
            continue

        generation_data = []
        with open(gen_file, "r") as f:
            for line in f:
                cmpd = Compound.from_tsv_line(line)
                generation_data.append(
                    {
                        "smiles": cmpd.smiles,
                        "id": cmpd.id,
                        "docking_score": cmpd.docking_score,
                    }
                )
        all_generations.append(generation_data)

    # Calculate global min/max scores and bin configuration
    all_scores = [
        compound["docking_score"]
        for generation in all_generations
        for compound in generation
    ]
    global_min = min(all_scores)
    global_max = max(all_scores)
    max_compounds_in_bin = 0
    bins = 15  # Fixed number of bins
    bin_width = (global_max - global_min) / bins

    # Calculate max count across all possible histogram configurations
    bin_edges = [global_min + i * bin_width for i in range(bins + 1)]
    for generation in all_generations:
        scores = [c["docking_score"] for c in generation]
        counts, _ = np.histogram(scores, bins=bin_edges)
        max_compounds_in_bin = max(max_compounds_in_bin, max(counts))

    # Add this line to get max generation size
    max_generation_size = min(
        50, max(len(generation) for generation in all_generations)
    )

    # Round up max_compounds_in_bin to nearest 5 or 10 for nice y-axis
    max_y = 5 * math.ceil(max_compounds_in_bin / 5)

    # Create the HTML template with embedded data
    html_template = f"""<!DOCTYPE html>
<html>
<head>
    <title>AutoGrow4 Results</title>
    <meta charset="UTF-8">
    <!-- Bootstrap CSS -->
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <!-- RDKit.js -->
    <script src="https://unpkg.com/smiles-drawer@2.0.1/dist/smiles-drawer.min.js"></script>
    <!-- Chart.js and Annotation plugin -->
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        .structure-container {{
            width: 150px;
            height: 150px;
            display: flex;
            align-items: center;
            justify-content: center;
        }}
        .structure-container canvas {{
            width: 150px !important;
            height: 150px !important;
        }}
        .slider-container {{
            margin: 20px 0;
        }}
    </style>
</head>
<body>
    <div class="container mt-4">
        <!-- Controls -->
        <div class="row slider-container">
            <div class="col-md-6">
                <label for="generationSlider" class="form-label">Generation: <span id="generationValue">1</span></label>
                <input type="range" class="form-range" id="generationSlider" min="1" max="{len(all_generations)}" value="1">
            </div>
            <div class="col-md-6">
                <label for="compoundsSlider" class="form-label">Top Compounds: <span id="compoundsValue">10</span></label>
                <input type="range" class="form-range" id="compoundsSlider" min="1" max="{max_generation_size}" value="10">
            </div>
        </div>

        <!-- Score Distribution Chart -->
        <div class="row mb-4">
            <div class="col">
                <div class="card">
                    <div class="card-body">
                        <!-- <h5 class="card-title">Score Distribution</h5> -->
                        <canvas id="histogram" style="max-height: 300px;"></canvas> <!-- Smaller height for above-the-fold view -->
                    </div>
                </div>
            </div>
        </div>

        <!-- Compounds Table -->
        <div class="row">
            <div class="col">
                <div class="card">
                    <div class="card-body">
                        <h5 class="card-title">Top Compounds</h5>
                        <div class="table-responsive">
                            <table class="table">
                                <thead>
                                    <tr>
                                        <th>Structure</th>
                                        <th>SMILES</th>
                                        <th>ID</th>
                                        <th>Docking Score</th>
                                    </tr>
                                </thead>
                                <tbody id="compoundsTable">
                                </tbody>
                            </table>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <script>
        // Embedded data
        const allGenerations = {json.dumps(all_generations)};
        const globalMin = {global_min};
        const globalMax = {global_max};
        const numBins = {bins};
        const binWidth = {bin_width};
        const maxY = {max_y};
        
        let smilesDrawer = new SmilesDrawer.Drawer({{
            width: 450,  // 3x larger for better resolution
            height: 450, // 3x larger for better resolution
            bondThickness: 2.0,  // Increase bond thickness for better visibility
            fontSizeLarge: 16,   // Increase font sizes proportionally
            fontSizeSmall: 14,
        }});

        let histogramChart = null;
        let moleculeQueue = [];
        let isProcessingQueue = false;

        // Initialize RDKit
        async function initRDKit() {{
            try {{
                // Wait for RDKit to be loaded
                if (!window.RDKit) {{
                    console.log("Waiting for RDKit...");
                    await new Promise(resolve => {{
                        window.addEventListener('RDKitReady', resolve);
                    }});
                }}
                console.log("RDKit loaded, initializing...");
                rdkit = await window.RDKit();
                console.log("RDKit initialized");
                processNextMolecule();
            }} catch (error) {{
                console.error('Failed to initialize RDKit:', error);
            }}
        }}
        
        window.addEventListener('DOMContentLoaded', function() {{
            initRDKit();  // Simplified initialization
            updateVisualization();
        }});
        
        window.addEventListener('DOMContentLoaded', function() {{
            // Check if RDKit is defined
            if (window.RDKit) {{
                if (typeof window.RDKit.then === 'function') {{
                    // If RDKit has a .then method (Promise-based), initialize with it
                    window.RDKit.then(initRDKit);
                }} else {{
                    // If RDKit is loaded directly without .then, initialize directly
                    initRDKit();
                }}
            }} else {{
                // If RDKit is not yet defined, listen for RDKitReady event
                window.addEventListener('RDKitReady', initRDKit);
            }}
            updateVisualization();
        }});

        async function generateMoleculeImage(smiles, containerId) {{
            if (!rdkit) {{
                console.log("RDKit not ready, queueing molecule:", smiles);
                moleculeQueue.push({{ smiles, containerId }});
                return;
            }}
            
            try {{
                console.log("Generating molecule for:", smiles);
                const mol = rdkit.get_mol(smiles);
                if (mol) {{
                    const svg = mol.get_svg();
                    const container = document.getElementById(containerId);
                    if (container) {{
                        container.innerHTML = svg;
                    }}
                    mol.delete();  // Clean up the molecule object
                }}
            }} catch (error) {{
                console.error('Error generating molecule:', error, smiles);
            }}
        }}

        async function processNextMolecule() {{
            if (!rdkit || isProcessingQueue || moleculeQueue.length === 0) {{
                if (!rdkit) console.log("RDKit not ready");
                if (isProcessingQueue) console.log("Already processing queue");
                if (moleculeQueue.length === 0) console.log("Queue empty");
                return;
            }}
            
            isProcessingQueue = true;
            console.log("Processing next molecule in queue");
            while (moleculeQueue.length > 0) {{
                const {{ smiles, containerId }} = moleculeQueue.shift();
                await generateMoleculeImage(smiles, containerId);
            }}
            isProcessingQueue = false;
        }}

        function calculateHistogramData(scores) {{
            const binCounts = Array(numBins).fill(0);
            const binLabels = [];
            
            for (let i = 0; i < numBins; i++) {{
                const binStart = globalMin + (i * binWidth);
                binLabels.push(binStart.toFixed(1));
            }}
            
            scores.forEach(score => {{
                const binIndex = Math.min(Math.floor((score - globalMin) / binWidth), numBins - 1);
                if (binIndex >= 0) {{
                    binCounts[binIndex]++;
                }}
            }});
            
            return {{ counts: binCounts, labels: binLabels }};
        }}

        function updateVisualization() {{
            const generationIndex = parseInt(document.getElementById('generationSlider').value);
            const topN = parseInt(document.getElementById('compoundsSlider').value);
            
            document.getElementById('generationValue').textContent = generationIndex;
            document.getElementById('compoundsValue').textContent = topN;
            
            const generationData = allGenerations[generationIndex - 1].slice(0, topN);

            const scores = generationData.map(c => c.docking_score);
            const avg = scores.reduce((a, b) => a + b) / scores.length;
            const median = scores.sort((a,b) => a-b)[Math.floor(scores.length/2)];

            console.log(avg, median);
            
            const histData = calculateHistogramData(scores);
            if (histogramChart) {{
                histogramChart.destroy();
            }}
            histogramChart = new Chart(document.getElementById('histogram'), {{
                type: 'bar',
                data: {{
                    labels: histData.labels,
                    datasets: [{{
                        label: 'Count',
                        data: histData.counts,
                        backgroundColor: '#8884d8'
                    }}]
                }},
                options: {{
                    responsive: true,
                    animation: false,  // Disable animation for immediate updates
                    plugins: {{
                        title: {{
                            display: true,
                            text: [`Mean: ${{avg.toFixed(2)}}  |  Median: ${{median.toFixed(2)}}`],
                            padding: {{
                                bottom: 10
                            }}
                        }},
                        legend: {{
                            display: false
                        }}
                    }},
                    scales: {{
                        y: {{
                            beginAtZero: true,
                            max: maxY,
                            title: {{
                                display: true,
                                text: 'Count'
                            }},
                            ticks: {{
                                stepSize: 1
                            }}
                        }},
                        x: {{
                            title: {{
                                display: true,
                                text: 'Docking Score'
                            }}
                        }}
                    }}
                }}
            }});
            
            const topCompounds = generationData
                .sort((a, b) => a.docking_score - b.docking_score)
                .slice(0, topN);

            const tableHtml = topCompounds.map((compound, index) => `
                <tr>
                    <td>
                        <div id="mol-${{index}}" class="structure-container">
                            <canvas width="150" height="150"></canvas>
                        </div>
                    </td>
                    <td class="text-truncate" style="max-width: 200px;">
                        ${{compound.smiles}}
                    </td>
                    <td class="text-truncate" style="max-width: 200px;">${{compound.id}}</td>
                    <td>${{compound.docking_score.toFixed(2)}}</td>
                </tr>
            `).join('');

            document.getElementById('compoundsTable').innerHTML = tableHtml;
            
            // Draw molecules
            topCompounds.forEach((compound, index) => {{
                const container = document.getElementById(`mol-${{index}}`);
                const canvas = container.querySelector('canvas');
                // Set high resolution canvas size
                canvas.width = 450;  // 3x larger
                canvas.height = 450; // 3x larger
                SmilesDrawer.parse(compound.smiles, function(tree) {{
                    smilesDrawer.draw(tree, canvas, 'light', false);
                }});
            }});
        }}

        document.getElementById('generationSlider').addEventListener('input', updateVisualization);
        document.getElementById('compoundsSlider').addEventListener('input', updateVisualization);
    </script>
</body>
</html>
"""

    # Write the HTML file
    output_file = os.path.join(output_dir, "summary.html")
    with open(output_file, "w") as f:
        f.write(html_template)

    log_info(f"Summary HTML file saved to: {output_file}")


def generate_summary_txt(output_dir: str) -> None:
    """Generate text-based summary of best compounds across all generations.
    
    Creates a ranked summary file containing the best compounds from all generations,
    along with a merged SDF file of their 3D structures.
    
    Args:
        output_dir (str): Path to AutoGrow output directory containing generation_X folders
    """
    import os
    import glob
    from autogrow.utils.logging import log_info, log_debug, log_warning

    # Read all generation data
    all_compounds = []
    generation_dirs = sorted(glob.glob(os.path.join(output_dir, "generation_*")))

    # Track stats for debugging
    total_lines = 0
    compounds_with_sdf = 0
    compounds_with_valid_sdf = 0

    for gen_dir in generation_dirs:
        gen_file = os.path.join(gen_dir, f"{os.path.basename(gen_dir)}_ranked.smi")
        if not os.path.exists(gen_file):
            continue

        log_debug(f"Processing generation directory: {gen_dir}")
        with open(gen_file, "r") as f:
            for line in f:
                total_lines += 1

                cmpd = Compound.from_tsv_line(line)

                # parts = line.strip().split("\t")

                # From your example file, the format is:
                # SMILES, ID, _, docking_score, diversity_score, sdf_path
                # if len(parts) >= 6:  # Make sure we have all needed fields
                compound = {
                    "smiles": cmpd.smiles,
                    "id": cmpd.id,
                    "docking_score": cmpd.docking_score,
                    "diversity_score": cmpd.diversity_score,
                    "sdf_path": cmpd.sdf_path,
                    "generation": os.path.basename(gen_dir),
                }

                # Verify SDF exists and is valid
                if os.path.exists(compound["sdf_path"]):
                    compounds_with_sdf += 1
                    if os.path.getsize(compound["sdf_path"]) > 0:
                        with open(compound["sdf_path"], "r") as test_f:
                            content = test_f.read()
                            if "$$$$" in content:  # Check if it's a valid SDF
                                compounds_with_valid_sdf += 1
                            else:
                                log_warning(
                                    f"Found SDF but it appears invalid: {compound['sdf_path']}"
                                )
                                compound["sdf_path"] = None
                    else:
                        log_warning(f"Found empty SDF file: {compound['sdf_path']}")
                        compound["sdf_path"] = None
                else:
                    log_warning(f"SDF file not found: {compound['sdf_path']}")
                    compound["sdf_path"] = None

                all_compounds.append(compound)

    # Log debug info
    log_debug(f"Total lines processed: {total_lines}")
    log_debug(f"Compounds with SDF path: {compounds_with_sdf}")
    log_debug(f"Compounds with valid SDF: {compounds_with_valid_sdf}")

    # Sort by docking score
    all_compounds.sort(key=lambda x: x["docking_score"])

    # Create summary files
    summary_tsv = os.path.join(output_dir, "summary_ranked.tsv")
    summary_sdf = os.path.join(output_dir, "summary_ranked.sdf")

    # Write ranked summary text file
    with open(summary_tsv, "w") as f:
        # Write header
        f.write(
            "Rank\tSMILES\tID\tDocking Score\tDiversity Score\tGeneration\tSDF Path\n"
        )

        # Write compound data
        for i, compound in enumerate(all_compounds, 1):
            f.write(
                f"{i}\t{compound['smiles']}\t{compound['id']}\t"
                f"{compound['docking_score']}\t{compound['diversity_score']}\t"
                f"{compound['generation']}\t"
                f"{compound['sdf_path'] or ''}\n"
            )

    # Create merged SDF file if SDF paths are available
    sdf_written = 0
    with open(summary_sdf, "w") as outfile:
        for rank, compound in enumerate(all_compounds, 1):
            if compound["sdf_path"] and os.path.exists(compound["sdf_path"]):
                try:
                    with open(compound["sdf_path"], "r") as infile:
                        content = infile.read()
                        if content.strip():  # Check if there's actual content
                            if "$$$$" in content:  # Verify it's a valid SDF
                                # Strip the last "$$$$" if it exists
                                content = content.strip()
                                if content.endswith("$$$$"):
                                    content = content[:-4].strip()

                                # Copy SDF content without the terminator
                                outfile.write(content)
                                outfile.write("\n")

                                # Add custom properties to SDF
                                outfile.write(f">  <Rank>\n{rank}\n\n")
                                outfile.write(
                                    f">  <Generation>\n{compound['generation']}\n\n"
                                )
                                outfile.write(f">  <ID>\n{compound['id']}\n\n")
                                outfile.write(
                                    f">  <Docking_Score>\n{compound['docking_score']}\n\n"
                                )
                                outfile.write(
                                    f">  <Diversity_Score>\n{compound['diversity_score']}\n\n"
                                )
                                outfile.write("$$$$\n")  # SDF separator
                                sdf_written += 1
                                log_debug(
                                    f"Successfully wrote structure {rank} from {compound['sdf_path']}"
                                )
                            else:
                                log_warning(
                                    f"Invalid SDF format in {compound['sdf_path']}"
                                )
                except Exception as e:
                    log_warning(f"Error processing {compound['sdf_path']}: {str(e)}")

    log_info("Created summary files:")
    log_info(f"  Ranked compounds: {summary_tsv}")
    log_info(f"  3D structures: {summary_sdf} ({sdf_written} structures)")

    # If no structures were written, something went wrong
    if sdf_written == 0:
        log_warning("No structures were written to the SDF file!")
        # List a few example compounds to help debugging
        for i, compound in enumerate(all_compounds[:5]):
            log_debug(f"Example compound {i+1}:")
            log_debug(f"  ID: {compound['id']}")
            log_debug(f"  SDF path: {compound['sdf_path']}")
            if compound["sdf_path"]:
                log_debug(f"  SDF exists: {os.path.exists(compound['sdf_path'])}")
                if os.path.exists(compound["sdf_path"]):
                    log_debug(
                        f"  SDF size: {os.path.getsize(compound['sdf_path'])} bytes"
                    )
