seurat2scanpy = function(seurat, ann.X = "SCT-data", ann.raw.X = "RNA-counts", h5ad_path = NULL) {
  
  if (is.null(h5ad_path)) {
    stop("必须提供输出的h5ad的绝对路径!")
  } else {
    tmp1 = stringr::str_split(h5ad_path, "/")[[1]]
    tmp2 = paste(tmp1[-length(tmp1)], collapse = "/")
    if (!dir.exists(tmp2)) dir.create(tmp2, recursive = TRUE)
  }
  
  if (is.null(ann.X)) {
    stop("必须提供ann.X!")
  }
  
  get_coord_valid = function(seurat, image = NULL) {
    spatial.df = GetTissueCoordinates(seurat, image = image)
    spatial.df = spatial.df %>% 
      rotate.axis(x = "imagerow", y = "imagecol", numBarcode = 96, angle = 90)
    spatial.df$imagerow = round(spatial.df$imagerow, 3)
    spatial.df$imagecol = round(spatial.df$imagecol, 3)
    
    p1 = ggplot(spatial.df, mapping = aes(x = imagerow, y = imagecol)) + geom_point() + theme(aspect.ratio = 1)
    p2 = SpatialDimPlot(seurat, group.by = "orig.ident", images = image, crop = F, pt.size.factor = 0.17) + theme(aspect.ratio = 1)
    print(wrap_plots(p1, p2))
    
    return(spatial.df)
  }
  
  .regularise_df = function(df, drop_single_values = TRUE) {
    if (ncol(df) == 0) df[["name"]] <- rownames(df)
    if (drop_single_values) {
      k_singular <- sapply(df, function(x) length(unique(x)) == 1)
      if (sum(k_singular) > 0) {
        warning(
          paste("Dropping single category variables:"),
          paste(colnames(df)[k_singular], collapse = ", ")
        )
      }
      df <- df[, !k_singular, drop = F]
      if (ncol(df) == 0) df[["name"]] <- rownames(df)
    }
    return(df)
  }
  
  ann.X.assay = stringr::str_split(ann.X, "-", simplify = T)[1, 1]
  ann.X.slots = stringr::str_split(ann.X, "-", simplify = T)[1, 2]
  ann.X.mat = Seurat::GetAssayData(object = seurat, assay = ann.X.assay, slot = ann.X.slots)
  
  # save metadata table:
  seurat$barcode = colnames(seurat)
  obs = .regularise_df(seurat@meta.data, drop_single_values = F)
  seurat@assays[[ann.X.assay]]@meta.features["highly_variable"] = rownames(seurat@assays[[ann.X.assay]]@meta.features) %in% VariableFeatures(seurat)
  if (sum(seurat@assays[[ann.X.assay]]@meta.features["highly_variable"]) == 0) {
    stop("seurat对象未设置高变基因!")
  }
  var = .regularise_df(Seurat::GetAssay(seurat, assay = ann.X.assay)@meta.features, drop_single_values = F)
  
  obsm = NULL
  reductions = names(seurat@reductions)
  if (length(reductions) > 0) {
    obsm = sapply(reductions, function(name) as.matrix(Seurat::Embeddings(seurat, reduction = name)), simplify = F)
    names(obsm) = paste0("X_", tolower(names(seurat@reductions)))
  }
  
  anndata = reticulate::import("anndata", convert = FALSE)
  
  if (length(seurat@images) > 0) {
    samples = str_sort(unique(seurat$orig.ident), numeric = T)
    spatial.df = Reduce(rbind, map(samples, ~flip.axis(get_coord_valid(seurat, image = make.names(.)), x = "imagerow", y = "imagecol", numBarcode = 96, horizontal = T)))
    spatial.df = spatial.df[colnames(seurat), ]
    
    obsm[["spatial"]] = as.matrix(spatial.df)
    
    imread = reticulate::import("matplotlib.image", convert = FALSE)
    library(reticulate)
    uns = dict("spatial" = dict())
    # write spatial image
    for (sample in str_sort(unique(seurat$orig.ident))) {
      temp_path = paste("/home/biogenger/Downloads/", paste0(sample, ".png"), sep = "/")
      if (sample %in% unique(seurat$orig.ident)) {
        sample = make.names(sample)
        png::writePNG(seurat@images[[sample]]@image, temp_path)
      } else {
        stop(paste(sample, "of orig.ident does not exist in images!"))
      }
      uns["spatial"][sample] = dict()
      if ("E12.5" %in% sample) {
        uns["spatial"][sample]["scalefactors"] = dict()
        uns["spatial"][sample]["scalefactors"]["spot_diameter_fullres"] = 1
        uns["spatial"][sample]["scalefactors"]["tissue_hires_scalef"] = 1
        uns["spatial"][sample]["scalefactors"]["fiducial_diameter_fullres"] = 184.32
        uns["spatial"][sample]["scalefactors"]["tissue_lowres_scalef"] = 1
      } else {
        uns["spatial"][sample]["scalefactors"] = dict()
        uns["spatial"][sample]["scalefactors"]["spot_diameter_fullres"] = 1
        uns["spatial"][sample]["scalefactors"]["tissue_hires_scalef"] = 1
        uns["spatial"][sample]["scalefactors"]["fiducial_diameter_fullres"] = 96
        uns["spatial"][sample]["scalefactors"]["tissue_lowres_scalef"] = 1
      }
      uns["spatial"][sample]["images"] = dict()
      uns["spatial"][sample]["images"]["hires"] = imread$imread(temp_path)
      system(paste0("rm -rf ", temp_path))
    }
    
    adata = anndata$AnnData(
      X = Matrix::t(ann.X.mat),
      obs = obs,
      var = var,
      obsm = obsm,
      uns = uns
    )
  } else {
    adata = anndata$AnnData(
      X = Matrix::t(ann.X.mat),
      obs = obs,
      var = var,
      obsm = obsm
    )
  }
  
  if (!is.null(ann.raw.X)) {
    ann.raw.X.assay = stringr::str_split(ann.raw.X, "-", simplify = T)[1, 1]
    ann.raw.X.slots = stringr::str_split(ann.raw.X, "-", simplify = T)[1, 2]
    ann.raw.X.mat = Seurat::GetAssayData(object = seurat, assay = ann.raw.X.assay, slot = ann.raw.X.slots)
    seurat@assays[[ann.raw.X.assay]]@meta.features["highly_variable"] = rownames(seurat@assays[[ann.raw.X.assay]]@meta.features) %in% VariableFeatures(seurat)
    var = .regularise_df(Seurat::GetAssay(seurat, assay = ann.raw.X.assay)@meta.features, drop_single_values = F)
    ann.raw = anndata$AnnData(
      X = Matrix::t(ann.raw.X.mat),
      var = var
    )
    adata$raw = ann.raw
  }
  
  adata$write(h5ad_path, compression = "gzip")
  
  return("转换成功!")
}


