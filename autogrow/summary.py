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
        with open(gen_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:  # Ensure we have at least SMILES, ID, and score
                    generation_data.append({
                        'smiles': parts[0],
                        'id': parts[1],
                        'docking_score': float(parts[4])
                    })
        all_generations.append(generation_data)

    # Calculate global min/max scores and bin configuration
    all_scores = [compound['docking_score'] 
                 for generation in all_generations 
                 for compound in generation]
    global_min = min(all_scores)
    global_max = max(all_scores)
    max_compounds_in_bin = 0
    bins = 15  # Fixed number of bins
    bin_width = (global_max - global_min) / bins

    # Calculate max count across all possible histogram configurations
    bin_edges = [global_min + i * bin_width for i in range(bins + 1)]
    for generation in all_generations:
        scores = [c['docking_score'] for c in generation]
        counts, _ = np.histogram(scores, bins=bin_edges)
        max_compounds_in_bin = max(max_compounds_in_bin, max(counts))
        
    # Add this line to get max generation size
    max_generation_size = min(50, max(len(generation) for generation in all_generations))

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
    <script src="https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js"></script>
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
        .structure-container svg {{
            max-width: 100%;
            max-height: 100%;
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
                                        <th>ID</th>
                                        <th>SMILES</th>
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
        
        let rdkit = null;
        let histogramChart = null;
        let moleculeQueue = [];
        let isProcessingQueue = false;

        // Initialize RDKit
        async function initRDKit() {{
            try {{
                rdkit = await window.RDKit;
                processNextMolecule();
            }} catch (error) {{
                console.error('Failed to initialize RDKit:', error);
            }}
        }}
        
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
                moleculeQueue.push({{ smiles, containerId }});
                return;
            }}
            
            try {{
                const mol = rdkit.get_mol(smiles);
                if (mol) {{
                    const svg = mol.get_svg();
                    const container = document.getElementById(containerId);
                    if (container) {{
                        container.innerHTML = svg;
                    }}
                }}
            }} catch (error) {{
                console.error('Error generating molecule:', error);
            }}
        }}

        async function processNextMolecule() {{
            if (!rdkit || isProcessingQueue || moleculeQueue.length === 0) return;
            
            isProcessingQueue = true;
            const {{ smiles, containerId }} = moleculeQueue.shift();
            await generateMoleculeImage(smiles, containerId);
            isProcessingQueue = false;
            processNextMolecule();
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
                            <div class="spinner-border text-primary" role="status">
                                <span class="visually-hidden">Loading...</span>
                            </div>
                        </div>
                    </td>
                    <td>${{compound.id}}</td>
                    <td class="text-truncate" style="max-width: 200px;">
                        ${{compound.smiles}}
                    </td>
                    <td>${{compound.docking_score.toFixed(2)}}</td>
                </tr>
            `).join('');

            document.getElementById('compoundsTable').innerHTML = tableHtml;
            
            topCompounds.forEach((compound, index) => {{
                generateMoleculeImage(compound.smiles, `mol-${{index}}`);
            }});
        }}

        document.getElementById('generationSlider').addEventListener('input', updateVisualization);
        document.getElementById('compoundsSlider').addEventListener('input', updateVisualization);
    </script>
</body>
</html>
"""

    # Write the HTML file
    output_file = os.path.join(output_dir, 'results_visualization.html')
    with open(output_file, 'w') as f:
        f.write(html_template)
