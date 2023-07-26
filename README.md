Wait for Document after published!

1. STvis (tool for interactive spatial cropping)
```
shiny_st(seurat_obj, assay = "SCT", slot = "data", image = "test", python_env = "~/miniconda3/envs/daily/bin/python", script = "~/script/filter_pixel_AI.py")
```
![image](https://github.com/EddieLv/STvis/assets/61786787/0a7e13cf-8ee4-44d6-9dbb-63c5150bce96)

2. ST Neighboring Chat (tool for spatial cell-cell communication based on neighboring method)
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

https://api.star-history.com/svg?repos=star-history/star-history
