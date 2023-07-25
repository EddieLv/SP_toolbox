Wait for Document after published!

1. STvis
```
shiny_st(seurat_obj, assay = "SCT", slot = "data", image = "test", python_env = "~/miniconda3/envs/daily/bin/python", script = "~/script/filter_pixel_AI.py")
```
![image](https://github.com/EddieLv/STvis/assets/61786787/0a7e13cf-8ee4-44d6-9dbb-63c5150bce96)

3. ST Neighboring Chat
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
```
```
result.df = run_st_neighbor_chat_genger(srat.traj, neighbor.res = neighbor.res, lr.df = lr_pair.df, assay = "SCT", slot = "data", label.col = "celltype.re", celltype_pair = c("AT1.pre", "AT1"), nperm = 2000)
result.df %>% dplyr::filter(pval < 0.05 & min.expr > 0.1) %>% arrange(pval)
#  ligand receptor celltype_pair ligand_celltype_pair_mean receptor_celltype_pair_mean mean_expr  min.expr    pval
#1  Vegfa     Nrp2   AT1.pre_AT1                0.20449837                  0.37398822 0.2892433 0.2129601 0.00375
#2   Fgf1    Fgfr2   AT1.pre_AT1                0.09344202                  1.76577291 0.9296075 0.1065280 0.01100
#3 Sema3a   Plxna1   AT1.pre_AT1                1.16676683                  0.09644369 0.6316053 0.1120997 0.02375
#4  Lamc2    Itga3   AT1.pre_AT1                0.33471254                  0.17916934 0.2569409 0.1286032 0.03000
```
```
srat.merge$Vegfa.expr = FetchData(srat.merge, vars = "Vegfa", slot = "data")[[1]]
srat.merge$Nrp2.expr = FetchData(srat.merge, vars = "Nrp2", slot = "data")[[1]]
Spatial2Featureplot_genger(srat.merge, orig.ident = "E16.5_slice18", image = "E16.5_slice18", show.image = F, featureA = "Vegfa.expr", featureB = "Nrp2.expr", show.label = T, theme.dark = F)
```
![image](https://github.com/EddieLv/STvis/assets/61786787/645b535a-0f5f-49eb-93fd-e82ad584dfe8)
