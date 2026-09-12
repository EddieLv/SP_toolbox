# ----------------------------
# Unpublished work!
# ----------------------------

# Installation
```
devtools::install_github("EddieLv/SP_toolbox")
```

# Introduction
Here I developed 3 functions for DBIT-seq data:
1.shiny_st() # Edit your spot in web freely!
2.seurat2scanpy() # Transform seurat.object to scanpy.anndata totally (including counts, data, dimreduc, neighbor, spaitla, etc.)!
3.run_st_neighbor_genger() # Cell communication algorithm based on spatial neighbor!

# 1. shiny_st (tool for interactive spatial cropping)

An interactive Shiny app for manually annotating, cropping, and labeling spatial transcriptomics spots with a tissue image as background.

## Function signature

```r
seurat_result <- shiny_st(
  seurat,                  # Seurat object (required)
  assay      = "SCT",      # Assay to use for gene expression query
  slot       = "data",     # Slot to use ("data", "counts", "scale.data")
  image      = NULL,       # Name of the image in seurat@images (required, e.g. "sample1")
  python_env = NULL,       # Path to Python executable for AI filtering (optional)
  script     = NULL,       # Path to filter_pixel_AI.py for AI filtering (optional)
  tooltip    = NULL,       # meta.data column name to show as tooltip on hover (optional)
  crop       = TRUE        # Crop the view to cells retained in the object
)
```

> **Return value**: a modified Seurat object. The app runs in the foreground — assign the result after quitting.

## Quick start

```r
library(SP.toolbox)

# Check the image name in your Seurat object
Images(seurat_obj)
# [1] "sample1"

# Launch the app and save the result
seurat_obj <- shiny_st(seurat_obj, assay = "SCT", slot = "data", image = "sample1")
```

For a slide containing several tissues, subset the cells first. With the default
`crop = TRUE`, `shiny_st()` frames only the retained tissue, matching the behavior
of `SpatialPlot(crop = TRUE)`:

```r
one_tissue <- subset(seurat_obj, cells = cells_from_one_tissue)
one_tissue <- shiny_st(one_tissue, image = "sample1", crop = TRUE)
```

Set `crop = FALSE` to keep showing the complete slide image.

## Step-by-step workflow

### Step 1 — Launch and verify the background image

After calling `shiny_st()`, a browser window opens showing the tissue H&E image with spots overlaid.
The sidebar on the left contains all controls.

### Step 2 — Color spots by feature or gene

| Control | Description |
|---|---|
| **Select feature** | Color spots by any `meta.data` column (e.g. `orig.ident`, `seurat_clusters`, `celltype`) |
| **Select shape** | Switch between square (`22`) and circle (`21`) spot shapes |
| **Gene input + √** | Type a gene name and click √ to color spots by its expression level |

### Step 3 — Align the tissue image (if spots are misaligned)

Use the image alignment controls in the sidebar to register spots onto the tissue:

| Control | Range | Description |
|---|---|---|
| **Flip by vertical** | — | Mirror spots top ↔ bottom |
| **Flip by horizontal** | — | Mirror spots left ↔ right |
| **Rotate by angle + √** | −360 to 360 | Rotate all spots around the center |
| **Move spots horizontally + √** | −96 to 96 | Shift all spots left / right |
| **Move spots vertically + √** | −96 to 96 | Shift all spots up / down |
| **Shrink spots horizontally + √** | 0 to 5 | Scale spot positions along the x-axis |
| **Shrink spots vertically + √** | 0 to 5 | Scale spot positions along the y-axis |

### Step 4 — Adjust spot appearance

| Control | Range | Description |
|---|---|---|
| **Spot.alpha** | 0–1 | Spot transparency (0 = invisible, 1 = fully opaque) |
| **Spot.size** | 0–1 | Spot rendering size |

### Step 5 — Subset spots (optional)

Use **Select feature to subset** + the text input to display only certain groups before lasso selection.
Type comma-separated values (e.g. `typeA,typeB`) and click **Select**.
Click **Back to all celltypes!** to restore all spots.

### Step 6 — Lasso select and label spots

1. Click the **lasso** icon in the top-right toolbar of the plot panel.
2. Draw a freehand boundary around the spots to annotate.
3. In **Set label for selected spots**, type a label (e.g. `tumor_region`).
4. Click **Confirm** — the label is written to the currently selected feature column.
5. Repeat Steps 1–4 for additional regions or labels.

> Tip: set **Select feature** to the column you are annotating so newly labeled spots change color immediately.

### Step 7 — (Optional) AI-based tissue filtering

Requires `python_env` and `script` to be provided at launch:

```r
seurat_obj <- shiny_st(
  seurat_obj,
  image      = "sample1",
  python_env = "~/miniconda3/envs/daily/bin/python",
  script     = "~/Biosoftwares/SP_toolbox/filter_pixel_AI.py"
)
```

In the app:
1. Adjust the **threshold** slider (0–100) to control filtering sensitivity.
2. Click **Power BY AI** — background spots are labeled `filtered`, tissue spots are labeled `exist` in a new column `ai.filter`.
3. Set **Filter mode** to `on` to switch spots to outline-only display for easier checking.

### Step 8 — Zoom and download

The toolbar in the top-right corner of the plot provides:

| Button | Function |
|---|---|
| **lasso** | Freehand spot selection mode |
| **zoom** | Zoom into a rectangular region |
| **recover** | Reset zoom to full view |
| **download png** | Save the current plot view as a PNG file |

### Step 9 — Quit and save

Click **Quit** in the sidebar. The app closes and returns the modified Seurat object with:

- All **Confirm**-ed labels written back into the corresponding `meta.data` column.
- An `ai.filter` column (`exist` / `filtered`) if AI filtering was used.

```r
# Inspect annotations after quitting
table(seurat_obj$orig.ident)
table(seurat_obj$ai.filter)   # only if AI filtering was run

# Save
saveRDS(seurat_obj, "seurat_annotated.rds")
```

## Full example

```r
library(SP.toolbox)

seurat_obj <- shiny_st(
  seurat_obj,
  assay      = "SCT",
  slot       = "data",
  image      = "sample1",
  python_env = "~/miniconda3/envs/daily/bin/python",
  script     = "~/Biosoftwares/SP_toolbox/filter_pixel_AI.py",
  tooltip    = "celltype"    # show celltype label on hover
)

saveRDS(seurat_obj, "seurat_annotated.rds")
```

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Background image not showing | ggplot2 ≥ 4.0 breaking change | Update: `devtools::install_github("EddieLv/SP_toolbox")` |
| Spots misaligned with image | Seurat ≥ 5.0 coordinate scaling change | Same — update the package |
| `Please set image!` error | `image` argument missing or wrong name | Run `Images(seurat_obj)` to get the exact name |
| AI button has no effect | Wrong `python_env` or `script` path | Verify both paths exist and the Python env has required packages |

![image](https://github.com/EddieLv/STvis/assets/61786787/0a7e13cf-8ee4-44d6-9dbb-63c5150bce96)

# 2. ST Neighboring Chat (tool for spatial cell-cell communication based on neighboring method)
> (repackaged from https://github.com/JohnGenome/ST-mouse-kidney-development)
```
sample = "E16.5_slice18"
srat.merge = readRDS("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/All/Paper/Fig4/srat.merge.reanno.rds")
srat.merge = SubsetSTData(srat.merge, spots = rownames(srat.merge@meta.data)[srat.merge$orig.ident == sample])
srat.merge$celltype.re = droplevels(srat.merge$celltype.re)
table(srat.merge$celltype.re)
#Sox9.pro   AT1.pre       AT1   AT2.pre       AT2 Matrix.FB    myo.FB      gCap 
#      55       802       103       690       226       955       576       646 
```
```
if (! dir.exists(ccc_path)) dir.create(ccc_path)
neighbor.res = run_st_neighbor_genger(srat.merge, coordx = "barcodeB", coordy = "barcodeA", label.col = "celltype.re", self.interaction = T, permut = T, nperm = 2000, max.distance = 2)
st_neighbor_heat_genger(neighbor.res$neighbor.df)
```
![image](https://github.com/EddieLv/STvis/assets/61786787/559946f8-a41d-47c6-8f31-09428d9b9c85)

```
st_celltype_pair_neighbor_genger(srat.traj, neighbor.res = neighbor.res, image = "E16.5_slice18", celltype_pair = c("AT1.pre", "AT1"), label.col = "celltype.re", pt.size = 2, show.image = T, show.label = T)
```
![image](https://github.com/EddieLv/STvis/assets/61786787/d2823afa-17d5-44e6-8c6a-db3cbc6e8c9c)

```
### Load LR-DB and start CCC Analysis  ###
lr_pair.df = read.table("/media/biogenger/D/LR_database/mouse_ligand_receptors_giotto.txt", header = T, stringsAsFactors = F)
lr_pair.df %>% head()
#  ligand receptor
#1    A2m     Lrp1
#2  Aanat   Mtnr1a
#3 Adam12    Itga9
#4 Adam12    Itgb1
#5 Adam12     Sdc4
#6 Adam15    Itga5
corrected_umi.ls = apply(srat.traj@assays$SCT@data, 1, sum)
corrected_umi.ls = corrected_umi.ls[corrected_umi.ls > 0]
lr_pair.df$mouseLigand.exist = ifelse(lr_pair.df$mouseLigand %in% names(corrected_umi.ls), "ligand", "no")
lr_pair.df$mouseReceptor.exist = ifelse(lr_pair.df$mouseReceptor %in% names(corrected_umi.ls), "receptor", "no")
lr_pair.df$lr.exist = paste0(lr_pair.df$mouseLigand.exist, "_", lr_pair.df$mouseReceptor.exist)
table(lr_pair.df$lr.exist)
lr_pair.df = lr_pair.df %>% dplyr::filter(lr.exist == "ligand_receptor")
lr_pair.df = lr_pair.df[ , 1:2]
colnames(lr_pair.df) = c("ligand", "receptor")
lr_pair.df$LR = paste0(lr_pair.df$ligand, "_", lr_pair.df$receptor)
```
```
result.df = run_st_neighbor_chat_genger(srat.merge, neighbor.res = neighbor.res, lr.df = lr_pair.df, assay = "SCT", slot = "data", label.col = "celltype.re", celltype_pair = c("AT1.pre", "myo.FB"), nperm = 2000)
result.df %>% head()
result.df %>% dplyr::filter(pval.both < 0.05 & min.expr.both > 0.1) %>% arrange(pval.both)
result.df %>% dplyr::filter(pval.ind < 0.05 & min.expr.ind > 0.1) %>% arrange(pval.ind)
#  ligand receptor           LR  celltype_pair ligand_celltype_both_mean receptor_celltype_both_mean ligand_celltype1_mean receptor_celltype2_mean min.expr.both min.expr.ind pval.ligand.both pval.receptor.both
#1    A2m     Lrp1     A2m_Lrp1 AT1.pre_myo.FB               0.022416981                  0.30527435           0.019948741               0.3688052    0.02241698   0.01994874           0.9950             0.9720
#2  Aanat   Mtnr1a Aanat_Mtnr1a AT1.pre_myo.FB               0.004674503                  0.00000000           0.006085798               0.0000000    0.00000000   0.00000000           0.3920             0.0000
#3 Adam12    Itga9 Adam12_Itga9 AT1.pre_myo.FB               0.821827418                  0.78610743           0.796703430               0.9556919    0.78610743   0.79670343           1.0000             0.0000
#4 Adam12    Itgb1 Adam12_Itgb1 AT1.pre_myo.FB               0.821827418                  0.82841966           0.796703430               0.8369752    0.82182742   0.79670343           1.0000             0.9985
#5 Adam12     Sdc4  Adam12_Sdc4 AT1.pre_myo.FB               0.821827418                  0.18802907           0.796703430               0.1787223    0.18802907   0.17872231           1.0000             0.6115
#6 Adam15    Itga5 Adam15_Itga5 AT1.pre_myo.FB               0.035670005                  0.08477808           0.033044533               0.1134542    0.03567000   0.03304453           0.8035             1.0000
#  pval.both pval.ligand.ind pval.receptor.ind pval.ind
#1   0.98350           0.989            0.0010  0.49500
#2   0.19600           0.154            0.0000  0.07700
#3   0.50000           1.000            0.0000  0.50000
#4   0.99925           1.000            0.8695  0.93475
#5   0.80575           1.000            0.7785  0.88925
#6   0.90175           0.854            0.3355  0.59475
```
```
srat.merge$Vegfa.expr = FetchData(srat.merge, vars = "Vegfa", slot = "data")[[1]]
srat.merge$Nrp2.expr = FetchData(srat.merge, vars = "Nrp2", slot = "data")[[1]]
Spatial2Featureplot_genger(srat.merge, orig.ident = "E16.5_slice18", image = "E16.5_slice18", show.image = F, featureA = "Vegfa.expr", featureB = "Nrp2.expr", show.label = T, theme.dark = F)
```
![image](https://github.com/EddieLv/STvis/assets/61786787/645b535a-0f5f-49eb-93fd-e82ad584dfe8)

# Stargazers
	
[![Stargazers over time](https://starchart.cc/EddieLv/STvis.svg)](https://starchart.cc/EddieLv/STvis)

<br><a href="https://github.com/Charmve/computer-vision-in-action#-以用促学先会后懂-"><img align="right" alt="Go for it!" src="https://raw.githubusercontent.com/Charmve/computer-vision-in-action/dd292873828228a753a9bd2de4576dbf8cc3902c/res/ui/footer-rocket.svg" height="220" title="Do what you like, and do it best!"/></a>
<br>
<p align="center">Feel free to ask any questions, open a PR if you feel something can be done differently!</p>
<h2 align="center">🌟 Star this repository 🌟</h2>
<p align="center">Created by <a href="https://github.com/EddieLv">EddieLv</a></p>
