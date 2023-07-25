### ST neighbor ###
get_neighbor_count = function(seurat, dist.df = NULL, label.col = NULL, permut = F, max.distance = max.distance) {
  # celltype permutation
  if(permut) {
    cells.celltype = sample(as.character(seurat@meta.data[[label.col]]))
  } else {
    cells.celltype = as.character(seurat@meta.data[[label.col]])
  }
  cells.celltype.dic = as.character(cells.celltype)
  names(cells.celltype.dic) = rownames(seurat@meta.data)

  dist.df.final = dist.df %>% mutate("cell1.type" = unname(cells.celltype.dic[cell1]), "cell2.type" = unname(cells.celltype.dic[cell2]))
  # max distance between two points
  dist.df.final = dist.df.final[dist.df.final$dist < max.distance, ]

  # output
  neighbor.df = as.data.frame(table(dist.df.final$cell1.type, dist.df.final$cell2.type))
  neighbor.df$Var1 = as.character(neighbor.df$Var1)
  neighbor.df$Var2 = as.character(neighbor.df$Var2)
  colnames(neighbor.df) = c("cell1.type", "cell2.type", "n.neighbor")

  return(list("neighbor.df" = neighbor.df, "dist.df.final" = dist.df.final))
}

get_celltype_neighbor = function(seurat, coordx = "barcodeB", coordy = "barcodeA", label.col = NULL, self.interaction = F, permut = F, nperm = 2000, max.distance = max.distance) {
  # get a symmetric matrix reflecting the euclidean distance between each two cells
  distx.m = seurat@meta.data[[coordx]] %*% matrix(1, nrow = 1, ncol = dim(seurat@meta.data)[1]) - matrix(1, nrow = dim(seurat@meta.data)[1], ncol = 1) %*% t(seurat@meta.data[[coordx]])
  disty.m = seurat@meta.data[[coordy]] %*% matrix(1, nrow = 1, ncol = dim(seurat@meta.data)[1]) - matrix(1, nrow = dim(seurat@meta.data)[1], ncol = 1) %*% t(seurat@meta.data[[coordy]])
  # euclidean distance
  dist.m = sqrt(distx.m**2 + disty.m**2)
  rm(list = c("distx.m", "disty.m"))
  # control self-self interaction
  if (!self.interaction) {
    dist.m = dist.m + diag(9999, nrow = dim(seurat@meta.data)[1])
  }
  rownames(dist.m) = rownames(seurat@meta.data)
  colnames(dist.m) = rownames(seurat@meta.data)
  # transform dist.m to long df
  dist.df = dist.m %>%
    reshape2::melt(varnames = c("cell1", "cell2"), value.name = "dist")
  dist.df$cell1 = as.character(dist.df$cell1)
  dist.df$cell2 = as.character(dist.df$cell2)
  message("Successfully create a cell-cell distance matrix!")

  # Get celltype neighboring info
  neighbor.res = get_neighbor_count(seurat, dist.df = dist.df, label.col = label.col, permut = F, max.distance = max.distance)
  neighbor.df = neighbor.res$neighbor.df
  dist.df.final = neighbor.res$dist.df.final
  if (permut) {
    message("Start permutation, it is very time consuming, 2000 permutations about 5 minutes!")
    random.m = mclapply(1:nperm, function(i) {get_neighbor_count(seurat, dist.df = dist.df, label.col = label.col, permut = T, max.distance = max.distance)$neighbor.df$n.neighbor}, mc.cores = 6, mc.preschedule = T)
    random.m = Reduce(cbind, random.m)
    random.m = as.matrix(random.m)
  } else {
    random.m = NULL
  }

  return(list("neighbor.df" = neighbor.df, "perm.mat" = random.m, "dist.df.final" = dist.df.final))
}

run_st_neighbor_genger = function(seurat, coordx = "barcodeB", coordy = "barcodeA", label.col = NULL, self.interaction = F, permut = T, nperm = 2000, max.distance = 2) {
  library(parallel)
  library(dplyr)
  library(reshape2)

  seurat@meta.data[[label.col]] = droplevels(seurat@meta.data[[label.col]])

  if (is.factor(seurat@meta.data[[label.col]])) {
    celltype.levels = levels(seurat@meta.data[[label.col]])
  } else {
    celltype.levels = unique(seurat@meta.data[[label.col]])
  }
  message(paste("celltypes within", paste(celltype.levels, collapse = ", "), "will be analyzed!"))

  celltype.neighbor.lis = get_celltype_neighbor(seurat, coordx = coordx, coordy = coordy, label.col = label.col, self.interaction = self.interaction, permut = permut, nperm = nperm, max.distance = max.distance)
  neighbor.df = celltype.neighbor.lis$neighbor.df
  random.m = celltype.neighbor.lis$perm.mat
  dist.df.final = celltype.neighbor.lis$dist.df.final

  neighbor.df$expect = apply(random.m, 1, mean)
  neighbor.df$oe = neighbor.df$n.neighbor / neighbor.df$expect
  # permutation test
  # H0: n.neighbor of celltypeA-celltypeB is not enriched
  # p(celltypeA-celltypeB) = sum(n.neighbor - random.n.neighbor < 0) / nperm
  # p is the proportion of mistake when H0 is rejected, that is, the acception ratio of H0
  neighbor.df$pvalue = apply((neighbor.df$n.neighbor - random.m) < 0, 1, sum) / nperm

  levels(neighbor.df$cell1.type) = celltype.levels
  levels(neighbor.df$cell2.type) = celltype.levels

  return(list("neighbor.df" = neighbor.df, "dist.df.final" = dist.df.final))
}
# neighbor.res = run_st_neighbor_genger(srat.traj, coordx = "barcodeB", coordy = "barcodeA", label.col = "celltype.re", self.interaction = T, permut = T, nperm = 2000, max.distance = 2)

st_neighbor_heat_genger = function(neighbor.df) {
  library(ggplot2)
  library(ggh4x)
  library(reshape2)

  neighbor.mat = neighbor.df%>% reshape2::acast(cell1.type ~ cell2.type, value.var = "n.neighbor")

  yclust = hclust(dist(neighbor.mat, method = "euclidean"), method = "ward.D2")
  xclust = hclust(dist(t(neighbor.mat), method = "euclidean"), method = "ward.D2")

  neighbor.df$oe.above1 = factor(ifelse(neighbor.df$oe > 1, "True", "False"), levels = c("True", "False"))
  p = ggplot(neighbor.df, aes(x = cell1.type, y = cell2.type)) +
    geom_tile(aes(fill = n.neighbor)) +
    geom_point(aes(color = oe.above1)) +
    geom_text(aes(label = format(pvalue, scientific = T)), vjust = 2) +
    scale_fill_gradient2(midpoint = median(neighbor.df$n.neighbor), low = "blue", mid = "white", high = "red") +
    scale_color_manual(values = list("False" = "transparent", "True" = "black")) +
    scale_y_dendrogram(hclust = yclust) +
    scale_x_dendrogram(hclust = xclust,position = "bottom") +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          plot.background = element_rect(fill = NA),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())

  return(p)
}
# st_neighbor_heat_genger(neighbor.res$neighbor.df)

### ST neighbor visualize ###
get_neighbor_cells = function(neighbor.result, celltype_pair = NULL) {
  dist.df = neighbor.result$dist.df.final
  dist.df = neighbor.result$dist.df.final %>% dplyr::filter(dist.df$cell1.type == celltype_pair[1] & dist.df$cell2.type == celltype_pair[2])
  celltype.pair.cell.list = list(unique(dist.df$cell1), unique(dist.df$cell2))
  names(celltype.pair.cell.list) = celltype_pair

  dist.df = neighbor.result$dist.df.final
  celltype.neighbor.cells = unlist(celltype.pair.cell.list)
  dist.df = neighbor.result$dist.df.final %>% dplyr::filter(dist.df$cell2 %in% celltype.neighbor.cells)
  all.neighbor.cells= unique(dist.df$cell1)

  return(list("celltype.pair.cell.list" = celltype.pair.cell.list, "all.neighbor.cells" = all.neighbor.cells))
}

st_celltype_pair_neighbor_genger = function(seurat, neighbor.res, image = NULL, celltype_pair = c("celltypeA", "celltyeB"), label.col = NULL, pt.size = 2, show.image = F, show.label = T) {
  message("Attention: coordinates transformation is based on 96 barcode!")
  numBarcode = 96

  seurat@meta.data[[label.col]] = droplevels(seurat@meta.data[[label.col]])
  meta.df = seurat@meta.data %>% rownames_to_column(var = "cell")

  celltype.pair.neighbor.res = get_neighbor_cells(neighbor.res, celltype_pair)
  all.neighbor.cells = celltype.pair.neighbor.res$all.neighbor.cells
  celltype.pair.cell.list = celltype.pair.neighbor.res$celltype.pair.cell.list

  meta.df$type = "other"
  meta.df$type[meta.df$cell %in% all.neighbor.cells & meta.df$celltype.re == celltype_pair[1]] = celltype_pair[1]
  meta.df$type[meta.df$cell %in% all.neighbor.cells & meta.df$celltype.re == celltype_pair[2]] = celltype_pair[2]

  meta.df$neighbor = "other"
  meta.df$neighbor[meta.df$cell %in% all.neighbor.cells] = "neighbor"

  type.cols = list("grey", "green", "blue")
  names(type.cols) = c("other", celltype_pair)

  neighbor.shapes = list(15, 17)
  names(neighbor.shapes) = c("other", "neighbor")

  img = seurat@images[[image]]@image
  img_grob = grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))

  coordinates = data.frame(iB = rep(seq(1, numBarcode, 1), each = length(seq(1, numBarcode, 1))),
                           iA = rep(seq(1, numBarcode, 1), times = length(seq(1, numBarcode, 1)))) %>%
    mutate(imagecol = (numBarcode + 1 - 0.5 - iA) * 1080 / numBarcode,
           imagerow = (iB - 0.5) * 1080 / numBarcode)

  coordinates$cell = paste0(coordinates$iB, "x", coordinates$iA)
  coordinates[ , c("imagerow", "imagecol")] = coordinates[ , c("imagerow", "imagecol")] %>%
    mutate(imagerow = imagerow * seurat@images[[image]]@scale.factors$lowres,
           imagecol = imagecol * seurat@images[[image]]@scale.factors$lowres) %>%
    rotate.axis.shiny(x = "imagerow", y = "imagecol", numBarcode = numBarcode, angle = 90) %>%
    rotate.axis.shiny(x = "imagerow", y = "imagecol", numBarcode = numBarcode, vertical = T)

  clean.barcodes = str_match(meta.df$cell, pattern = "\\d+x\\d+")[ , 1]
  prefix = gsub(clean.barcodes[1], "", meta.df$cell[1])
  rownames(meta.df) = gsub(prefix, "", meta.df$cell)
  meta.df[coordinates$cell, c("X", "Y")] = coordinates[ , c("imagerow", "imagecol")]

  if (show.image) {
    p = ggplot(meta.df, aes(x = X, y = Y)) +
      annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
      geom_point(aes(color = type, shape = neighbor), size = pt.size) +
      scale_color_manual(values = type.cols) +
      scale_shape_manual(values = neighbor.shapes) +
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = T, clip = "on") +
      ylim(nrow(img), 0) + xlim(0, ncol(img)) +
      theme_void() +
      theme(legend.position = "right",
            aspect.ratio = 1)
  } else {
    meta.df[ , c("X", "Y")] = meta.df[ , c("X", "Y")] %>% flip.axis.shiny(x = "X", y = "Y", numBarcode = numBarcode, horizontal = T)
    p = ggplot(meta.df, aes(x = X, y = Y)) +
      geom_point(aes(color = type, shape = neighbor), size = pt.size) +
      scale_color_manual(values = type.cols) +
      scale_shape_manual(values = neighbor.shapes) +
      scale_x_continuous(expand = c(0, 0), breaks = sort(unique(round(coordinates$imagerow, 3))), labels = 1:numBarcode, limits = c(0, 1080), sec.axis = dup_axis()) +
      scale_y_continuous(expand = c(0, 0), breaks = sort(unique(round(coordinates$imagecol, 3))), labels = 1:numBarcode, limits = c(0, 1080), sec.axis = dup_axis()) +
      theme_void() +
      theme(legend.position = "right",
            aspect.ratio = 1,
            panel.border = element_blank(),
            axis.line.x = element_line(color = "black", linewidth = 0.5, linetype = "solid"),
            axis.line.y = element_line(color = "black", linewidth = 0.5, linetype = "solid"),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank())
  }

  if (show.label & !show.image) {
    p = p +
      theme(axis.ticks.length = unit(0.1, "cm"),
            axis.ticks.x = element_line(color = "black", linewidth = 0.3, linetype = "solid"),
            axis.ticks.y = element_line(color = "black", linewidth = 0.3, linetype = "solid"),
            axis.text.x = element_text(size = 8, colour = "black", angle = 270, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5, vjust = 0.5),
            axis.text.x.top = element_blank(),
            axis.ticks.x.top = element_blank(),
            axis.text.y.right = element_blank(),
            axis.ticks.y.right = element_blank())
  }

  return(p)
}
# st_celltype_pair_neighbor_genger(srat.traj, neighbor.res = neighbor.res, image = "E16.5_slice18", celltype_pair = c("AT1.pre", "AT1"), label.col = "celltype.re", pt.size = 2, show.image = T, show.label = T)

### ST neighbor CCC ###
celltype_permutate = function(meta.df, selected.cells = NULL, label.col = NULL, celltype_pair = NULL, permut = F) {
  meta.df = meta.df %>% dplyr::filter(cell %in% selected.cells)
  if (permut) {
    meta.df[[label.col]] = sample(meta.df[[label.col]], dim(meta.df)[1])
  }

  permut.res = split(meta.df$cell, meta.df[[label.col]])
  permut.res = permut.res[celltype_pair]

  return(permut.res)
}

get_LR_expr_score = function(count.m, nei_cells = NULL, lr.df = NULL) {
  lr.df$celltype_pair = paste(names(nei_cells), collapse = "_")

  Ligand_both_mean = apply(count.m[lr.df$ligand, unname(unlist(nei_cells))], 1, mean)
  Receptor_both_mean = apply(count.m[lr.df$receptor, unname(unlist(nei_cells))], 1, mean)
  lr.df$ligand_celltype_both_mean = Ligand_both_mean
  lr.df$receptor_celltype_both_mean = Receptor_both_mean

  Ligand_mean = apply(count.m[lr.df$ligand, nei_cells[[1]]], 1, mean)
  Receptor_mean = apply(count.m[lr.df$receptor, nei_cells[[2]]], 1, mean)
  lr.df$ligand_celltype1_mean = Ligand_mean
  lr.df$receptor_celltype2_mean = Receptor_mean

  return(lr.df)
}

run_st_neighbor_chat_genger = function(seurat, neighbor.res = NULL, lr.df = NULL, assay = "SCT", slot = "data", label.col = NULL, celltype_pair = c("celltypeA", "celltyeB"), nperm = 2000, both.expr = T) {
  DefaultAssay(seurat) = assay
  meta.df = seurat@meta.data %>% rownames_to_column(var = "cell")
  expr.mat = as.matrix(GetAssayData(seurat, slot = slot))

  celltype.pair.neighbor.res = get_neighbor_cells(neighbor.res, celltype_pair)
  all.neighbor.cells = celltype.pair.neighbor.res$all.neighbor.cells

  neighbor_LR_score = get_LR_expr_score(count.m = expr.mat,
                                        nei_cells = celltype_permutate(meta.df, selected.cells = all.neighbor.cells, label.col = label.col, celltype_pair = celltype_pair, permut = F),
                                        lr.df = lr_pair.df)
  message("Start permutation, 2000 permutations about 1 minute!")
  # consider LR expression on both celltypes
  neighbor_LR_score.permut = Reduce(function(i, j) {return(list(cbind(i[[1]], j[[1]]), cbind(i[[2]], j[[2]])))},
                                    mclapply(1:nperm, function(i) {
                                      neighbor_LR_score = get_LR_expr_score(count.m = expr.mat,
                                                                            nei_cells = celltype_permutate(meta.df, selected.cells = all.neighbor.cells, label.col = "celltype.re", celltype_pair = celltype_pair, permut = T),
                                                                            lr.df = lr_pair.df)
                                      res.list = list(neighbor_LR_score$ligand_celltype_both_mean,
                                                      neighbor_LR_score$receptor_celltype_both_mean)
                                      return(res.list)},
                                      mc.cores = 6))
  neighbor_lr_pv_both = cbind(apply((neighbor_LR_score.permut[[1]] - neighbor_LR_score$ligand_celltype_both_mean) > 0, 1, sum) / nperm,
                              apply((neighbor_LR_score.permut[[2]] - neighbor_LR_score$receptor_celltype_both_mean) > 0, 1, sum) / nperm)
  # consider LR expression on individual celltypes
  neighbor_LR_score.permut = Reduce(function(i, j) {return(list(cbind(i[[1]], j[[1]]), cbind(i[[2]], j[[2]])))},
                                    mclapply(1:nperm, function(i) {
                                      neighbor_LR_score = get_LR_expr_score(count.m = expr.mat,
                                                                            nei_cells = celltype_permutate(meta.df, selected.cells = all.neighbor.cells, label.col = "celltype.re", celltype_pair = celltype_pair, permut = T),
                                                                            lr.df = lr_pair.df)
                                      res.list = list(neighbor_LR_score$ligand_celltype1_mean,
                                                      neighbor_LR_score$receptor_celltype2_mean)
                                      return(res.list)},
                                      mc.cores = 6))
  neighbor_lr_pv_ind = cbind(apply((neighbor_LR_score.permut[[1]] - neighbor_LR_score$ligand_celltype1_mean) > 0, 1, sum) / nperm,
                             apply((neighbor_LR_score.permut[[2]] - neighbor_LR_score$receptor_celltype2_mean) > 0, 1, sum) / nperm)

  # add statistics
  neighbor_LR_score$min.expr.both = pmin(neighbor_LR_score$ligand_celltype_both_mean, neighbor_LR_score$receptor_celltype_both_mean)
  neighbor_LR_score$min.expr.ind = pmin(neighbor_LR_score$ligand_celltype1_mean, neighbor_LR_score$receptor_celltype2_mean)

  neighbor_lr_pv_both_mean = apply(neighbor_lr_pv_both, 1, mean)
  neighbor_LR_score[ , c("pval.ligand.both", "pval.receptor.both")] = neighbor_lr_pv_both
  neighbor_LR_score$pval.both = neighbor_lr_pv_both_mean

  neighbor_lr_pv_ind_mean = apply(neighbor_lr_pv_ind, 1, mean)
  neighbor_LR_score[ , c("pval.ligand.ind", "pval.receptor.ind")] = neighbor_lr_pv_ind
  neighbor_LR_score$pval.ind = neighbor_lr_pv_ind_mean

  return(neighbor_LR_score)
}
# result.df = run_st_neighbor_chat_genger(srat.traj, neighbor.res = neighbor.res, lr.df = lr_pair.df, assay = "SCT", slot = "data", label.col = "celltype.re", celltype_pair = c("AT1.pre", "AT1"), nperm = 2000)

# lr_pair.df = read.table("/media/biogenger/D/LR_database/mouse_ligand_receptors_giotto.txt", header = T, stringsAsFactors = F)
# corrected_umi.ls = apply(srat.traj@assays$SCT@data, 1, sum)
# corrected_umi.ls = corrected_umi.ls[corrected_umi.ls > 0]
# lr_pair.df$mouseLigand.exist = ifelse(lr_pair.df$mouseLigand %in% names(corrected_umi.ls), "ligand", "no")
# lr_pair.df$mouseReceptor.exist = ifelse(lr_pair.df$mouseReceptor %in% names(corrected_umi.ls), "receptor", "no")
# lr_pair.df$lr.exist = paste0(lr_pair.df$mouseLigand.exist, "_", lr_pair.df$mouseReceptor.exist)
# table(lr_pair.df$lr.exist)
# lr_pair.df = lr_pair.df %>% dplyr::filter(lr.exist == "ligand_receptor")
# lr_pair.df = lr_pair.df[ , 1:2]
# colnames(lr_pair.df) = c("ligand", "receptor")

