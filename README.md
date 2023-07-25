Wait for Document after published!

1. ST Neighboring Chat
> (repackaged from https://github.com/JohnGenome/ST-mouse-kidney-development)
```
sample = "E16.5_slice18"
srat.merge = readRDS("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/All/Paper/Fig4/srat.merge.reanno.rds")
srat.merge = SubsetSTData(srat.merge, spots = rownames(srat.merge@meta.data)[srat.merge$orig.ident == sample])
srat.merge$celltype.re = droplevels(srat.merge$celltype.re)
table(srat.merge$celltype.re)
#Sox9.pro   AT1.pre       AT1   AT2.pre       AT2 Matrix.FB    myo.FB      gCap 
#      55       802       103       690       226       955       576       646 
       
if (! dir.exists(ccc_path)) dir.create(ccc_path)
neighbor.res = run_st_neighbor_genger(srat.merge, coordx = "barcodeB", coordy = "barcodeA", label.col = "celltype.re", self.interaction = T, permut = T, nperm = 2000, max.distance = 2)
st_neighbor_heat_genger(neighbor.res$neighbor.df)
st_celltype_pair_neighbor_genger(srat.traj, neighbor.res = neighbor.res, image = "E16.5_slice18", celltype_pair = c("AT1.pre", "AT1"), label.col = "celltype.re", pt.size = 2, show.image = T, show.label = T)

### Load LR-DB and start CCC Analysis  ###
lr_pair.df = read.table("/media/biogenger/D/LR_database/mouse_ligand_receptors_giotto.txt", header = T, stringsAsFactors = F)
corrected_umi.ls = apply(srat.traj@assays$SCT@data, 1, sum)
corrected_umi.ls = corrected_umi.ls[corrected_umi.ls > 0]
lr_pair.df$mouseLigand.exist = ifelse(lr_pair.df$mouseLigand %in% names(corrected_umi.ls), "ligand", "no")
lr_pair.df$mouseReceptor.exist = ifelse(lr_pair.df$mouseReceptor %in% names(corrected_umi.ls), "receptor", "no")
lr_pair.df$lr.exist = paste0(lr_pair.df$mouseLigand.exist, "_", lr_pair.df$mouseReceptor.exist)
table(lr_pair.df$lr.exist)
lr_pair.df = lr_pair.df %>% dplyr::filter(lr.exist == "ligand_receptor")
lr_pair.df = lr_pair.df[ , 1:2]
colnames(lr_pair.df) = c("ligand", "receptor")

result.df = run_st_neighbor_chat_genger(srat.traj, neighbor.res = neighbor.res, lr.df = lr_pair.df, assay = "SCT", slot = "data", label.col = "celltype.re", celltype_pair = c("AT1.pre", "AT1"), nperm = 2000)
result.df %>% dplyr::filter(pval < 0.05 & min.expr > 0.1) %>% arrange(pval)

```
