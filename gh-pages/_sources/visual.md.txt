# Visualization

After performing the analyses, ***DiNetxify*** offers visualization tools to help interpret the results. The visualizations tie together the findings from PheWAS, comorbidity network, and disease trajectory analyses in a coherent way. The primary class for visualization is `dnt.visualization.Plot`, which takes the result DataFrames and generates interactive plots.

## Initializing the Plot object

To create visualizations, first initialize a `Plot` object with your analysis results. You will pass in the PheWAS, comorbidity network, and trajectory results DataFrames, and optionally specify how the “exposure” (the primary risk factor or cohort definition) should appear in the plots. For example:

```python
from DiNetxify.visualization import Plot  

# Suppose phewas_result, com_network_result, trajectory_result are obtained from above  
result_plot = Plot(  
    phewas_result=phewas_result,  
    comorbidity_result=com_network_result,  
    trajectory_result=trajectory_result,  
    exposure_name='Short LTL',    # Name of the exposure (for labeling the exposure node). Use None for exposed-only cohorts.  
    exposure_size=15,            # Relative size scaling for the exposure node (to make it prominent). None for exposed-only cohorts.  
    exposure_location=(0, 0, 0)  # 3D coordinates for the exposure node. If None, it defaults to (0,0,0).  
)  

# If this were an exposed-only cohort (no explicit exposure variable), you would set:  
result_plot = Plot(  
    phewas_result=phewas_result,  
    comorbidity_result=com_network_result,  
    trajectory_result=trajectory_result,  
    exposure_name=None,     # No separate exposure node  
    exposure_size=None,  
    exposure_location=None  
)  
```

- **phewas_result** – DataFrame from the PheWAS analysis (must include columns for phecode identifier, disease name/system, case counts, significance etc.).
- **comorbidity_result** – DataFrame from the comorbidity network analysis. It should include columns for disease pairs (D1, D2 identifiers), some unique pair name or ID, the association metrics (like beta or correlation), and a significance indicator (True/False for whether that pair is significant in the network).
- **trajectory_result** – DataFrame from the disease trajectory analysis. It should have columns for the disease pairs (D1, D2, typically oriented as source -> target), the effect sizes (like OR or HR), and a significance indicator for temporal association (True/False for adjusted p-value significance).
- **exposure_name** – A label for the exposure of interest (the factor that defines exposed vs unexposed in the cohort). In our examples, the exposure was “Short LTL” (short leukocyte telomere length) in one of the case studies, hence the example. If you are analyzing something like “smoking” or “diabetes” as the exposure, you’d put that. For an exposed-only study, use None (because there isn’t a separate exposure node to highlight).
- **exposure_location** – The (x, y, z) coordinates where the exposure node should be placed in the 3D plot. By default, if None, the exposure node will be placed at the origin (0,0,0). This is relevant only for 3D plotting; if exposure is None, this is ignored.
- **exposure_size** – A scaling factor for the exposure node’s size in the network visualization. Increase this to make the exposure node larger relative to disease nodes (to emphasize it). If None, in an exposed-only design, the exposure node is not present.

The `Plot` class will internally verify that the required columns exist in the input DataFrames (for example, it expects certain default column names like `'phecode_d1'`, `'phecode_d2'` for pair identifiers, `'phecode'` for disease codes in PheWAS, `'system'` for disease system category, `'name_disease_pair'` for a unique pair name, `'..._significance'` for significance flags, etc. If you did not change column names, it uses these defaults). You can override these defaults by passing optional parameters if your DataFrame uses different column names (we’ll mention those as needed in each plot function below).

Now that `result_plot` is created, we can generate specific plots. The `Plot` class provides multiple methods for different visualizations. Each produces either an interactive HTML file or a static image, as noted.

## PheWAS plot

The `result_plot.phewas_plot()` function creates a summary plot of the PheWAS results – essentially a visual depiction of which diseases were associated with the exposure. In a standard or matched cohort study, this is typically a Manhattan-style plot or bar plot showing hazard ratios (HRs) for each significant disease. In an exposed-only cohort, since we don’t have HRs, the plot might show the number of cases of each disease (to highlight which are most frequent). The plot also differentiates diseases by their category/system (often using color coding).

For example:

```python
# Generate a PheWAS plot  
result_plot.phewas_plot(  
    path="/your/project/path/phewas_plot.png",  # output file path (supports .png, .svg, .jpg)  
    is_exposure_only=False                     # False for cohort/matched designs; True if this is an exposed-only cohort  
)  
```

This will save an image to the specified path. You can use `.png` for a static image (good for publications), or `.svg` for a vector graphic, etc. If you want to just display it in a Jupyter notebook, you could omit the path and it might display inline, but typically you provide a path to save.

Parameters for `phewas_plot()` include:

- **path** – File path including filename and extension where the plot will be saved. Ensure the extension is one of the supported image types.
- **is_exposure_only** – Boolean flag; set to `True` if your analysis is an exposed-only cohort (so the plot knows it shouldn’t expect HRs and will plot counts instead). For our example (standard/matched cohort), it’s `False`.

**Optional parameters:** (if your DataFrame columns differ from defaults or if you want to adjust aesthetics)

- **col_coef** – Name of the column in `phewas_result` that contains the effect size (e.g., HR or OR). *(Default: 'phewas_coef')*.
- **col_se** – Name of the column for standard error of the effect size (used to plot error bars). *(Default: 'phewas_se')*.
- **col_system** – Column name for the disease system/category. *(Default: 'system')*.
- **col_disease** – Column name for the disease description/name. *(Default: 'disease')*.
- **col_exposure** – Column name for number of cases in exposed group (used in exposed-only plot). *(Default: 'N_cases_exposed')*.
- **disease_font_size** (parameter is actually `disese_font_size` due to a minor naming typo in code) – Font size for disease labels in the plot (if labels are shown). *(Default: 10)*.
- **system_font_size** – Font size for the disease system labels on the plot. *(Default: 17)*.
- **dpi** – Resolution of the output image (dots per inch). Higher DPI gives a higher resolution image. *(Default: 200)*.

Using these, the PheWAS plot will highlight which diseases came out as significantly associated. Typically, you’ll see something like a scatter of points or bars, colored by system, maybe labeled for the top hits. In our dummy data, since everything is random, the plot isn’t meaningful medically, but with real data this can quickly show you the pattern of associations.

## Comorbidity network plot

The `result_plot.comorbidity_network_plot()` function creates an interactive network visualization of the comorbidity relationships. It clusters diseases into communities based on their network connections (using the Louvain community detection algorithm) and plots them, often in a circular layout where each community occupies a sector. The nodes (diseases) are colored by their disease system, and edges represent significant associations from the comorbidity network analysis. This plot is typically output as an HTML file because it’s interactive (you can hover to see details, zoom, etc.).

For example:

```python
# Generate an interactive comorbidity network plot  
result_plot.comorbidity_network_plot(  
    path="/your/project/path/comorbidity_network.html"  
)  
```

This will save an HTML file which you can open in a web browser to explore. Each node (disease) might be labeled or have tooltip info (like disease name, maybe prevalence), and edges might have tooltips for correlation values. Nodes within the same community are grouped together. This visualization helps identify clusters of diseases that frequently co-occur beyond what would be expected.

Optional parameters for `comorbidity_network_plot()` include layout and styling options:

- **max_radius** – Maximum radial distance (in pixels) from the center that nodes can be placed. This essentially controls the size of the outer circle of the plot. *(Default: 90.0)*.
- **min_radius** – Minimum radial distance (pixels) from center for node placement (the inner boundary of the network). *(Default: 35.0)*.
- **layer_distance** – Radial distance between concentric layers of nodes (if nodes are layered by some criteria, e.g., significance or degree). *(Default: 40.0)*.
- **size_reduction** – A scaling factor (0 to 1) for node sizes to ensure they fit nicely (smaller values make nodes proportionally smaller). *(Default: 0.5)*.
- **line_width** – Width (pixels) of the lines (edges) connecting nodes. *(Default: 1.0)*.
- **line_color** – Color of the edges. Can specify a named color (e.g., `'steelblue'`), hex code (`'#4682B4'`), or RGB tuple. *(Default: 'black')*.

These parameters allow fine-tuning the appearance if needed (for example, if node labels overlap, you might reduce node sizes or adjust radii). Usually, the defaults produce a clear separation of communities.

## Disease trajectory plot

The `result_plot.trajectory_plot()` function generates visualizations for disease trajectories. Typically, it will produce one plot per community of diseases (as identified in the network analysis) to show how diseases progress within that community. For each significant disease pair (temporal relationship), an arrow or directed edge is drawn from the antecedent disease to the consequent disease. The output is often a set of static image files (e.g., one PNG per community) because it might be easier to print or inspect individually.

For example:

```python
# Generate disease trajectory plots (one per community)  
result_plot.trajectory_plot(  
    path="/your/project/path/trajectory_plots/"  
)  
```

> **Important:** Here, `path` is a directory (you should include the trailing slash). The function will then save multiple files in this directory, named perhaps by community or numbered sequentially. Ensure the directory exists or the function might attempt to create it.

Optional parameters for `trajectory_plot()` include:

- **source** – Column name in the DataFrames for the source disease (D1). *(Default: 'phecode_d1')*.
- **target** – Column name for the target disease (D2). *(Default: 'phecode_d2')*.
- **dpi** – Image resolution (like in PheWAS plot). *(Default: 500 for these plots, to ensure high clarity since many arrows/labels might be present.)*
- **cluster_weight** – This parameter specifies which edge weight from the network to use when arranging the layout. *(Default: 'comorbidity_beta')*, meaning it might use the beta coefficient from the comorbidity network as the weight for clustering layout. Usually, you won’t need to change this.

The trajectory plots will illustrate sequences. Within each community, you might see something like a directed acyclic graph (hopefully acyclic if the data suggests a direction). The arrow from disease A to B indicates A tends to precede B. By generating one plot per community, it keeps the graphs manageable and interpretable. If a community has no significant trajectories, it might be empty or skipped.

After running this, check the output directory for files. For instance, you might find files like `community_1.png`, `community_2.png`, etc., each showing the trajectory network for that community of diseases.

> *(Note: In interactive use, the function might print or log which communities are being plotted and the number of images saved.)*

## Three-dimensional plot

The `result_plot.three_dimension_plot()` function is a special visualization that combines the comorbidity network and trajectory information into a single 3D interactive figure. It essentially places the comorbidity network on one plane (say, the horizontal plane) and the trajectory connections on a vertical plane, giving a three-dimensional view where you can see both types of relationships simultaneously. The exposure node (if any) is usually at the center, and disease nodes are arranged around. This plot is typically an interactive HTML as well, since you’d want to rotate and zoom in 3D.

Example usage:

```python
# Generate a combined 3D network plot  
result_plot.three_dimension_plot(  
    path="/your/project/path/combined_network.html"  
)  
```

This saves an HTML file with the interactive 3D visualization. When you open it, you can drag to rotate the network in 3D space. One plane (e.g., viewed from top-down) will show the clusters and connections akin to the comorbidity network, and the other plane (viewed from the side) will show the directed trajectories. Where they intersect, you get a full picture of how diseases group and progress.

Optional parameters for `three_dimension_plot()` include various layout settings:

- **max_radius** – Maximum radius for node placement from the center (similar to the 2D network plot). *(Default: 180.0)* (since in 3D we might allow a larger radius).
- **min_radius** – Minimum radius from center. *(Default: 35.0)*.
- **layer_distance** – Distance between layers in the radial direction (similar concept to above). *(Default: 40.0)*.
- **layout_width** – Width of the overall figure in pixels. *(Default: 900.0)*.
- **layout_height** – Height of the figure in pixels. *(Default: 900.0)*.
- **line_color** – Color for trajectory lines (since in 3D plot, maybe comorbidity edges are one style and trajectory edges another). *(Default: 'black')*.
- **line_width** – Width of trajectory lines. *(Default: 1.0)*.
- **size_reduction** – Node size scaling factor (0.1–1.0). *(Default: 0.5)*.
- **cluster_reduction_ratio** – A factor (0.1–1.0) to compress or spread out clusters in the 3D space. Lower means clusters are more tightly grouped. *(Default: 0.4)*.
- **cluster_weight** – Which weight to use for determining cluster layout (as in trajectory_plot, usually 'comorbidity_beta'). *(Default: 'comorbidity_beta')*.
- **font_style** – Font family for text elements (node labels, etc.). *(Default: 'Times New Roman')*.
- **font_size** – Base font size for text. *(Default: 15.0)*.

The 3D plot is a bit advanced and may require a good computer/browser to manipulate if there are many nodes, but it provides a unique integrated view.

With all these visualization tools, even a beginner user can not only run the analysis with ***DiNetxify*** but also see the results in intuitive forms, which can greatly aid interpretation and presentation. Next, we provide a quick reference to the API of the main classes and functions for convenience.