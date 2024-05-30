create_seurat_RNA = function(count_mat = NULL, gene_id_input = F, metadata_path = NULL, 
                             assay = "RNA", project = "genger", species = "mouse", 
                             run_cellQC = F, min.nCount = 1000, max.nCount.quantile = 0.99, min.features = 500, max.nFeature.quantile = 0.99,
                             run_geneQC = T, min.cells = 5, genes_remove = NULL, 
                             cellname_prefix = NULL, cellname_suffix = NULL, 
                             isST = F, numBC = 96, image_dir = NULL, barcode_df = NULL, 
                             run_doublet = F) {
  ### 1.检查所有文件的合法性 ###
  if (isST) {
    # image_dir合法性
    if (!is.null(image_dir)) {
      if (!dir.exists(image_dir)) {
        stop(paste0("请检查image_dir:", image_dir))
      }
    }
    # barcode_df合法性
    if (!file.exists(barcode_df)) {
      stop(paste0("请检查barcode_df:", barcode_df))
    }
  }
  
  # count_mat合法性
  print("-------------- run --------------")
  message(paste0("输入的count_mat有", ncol(count_mat), "个细胞, ", nrow(count_mat), "个基因"))

  # metadata合法性
  if (!is.null(metadata_path)) {
    if (!file.exists(metadata_path)) {
      stop(paste0("请检查metadata_path:", metadata_path))
    } else {
      meta.data = fread(metadata_path, check.names = F)
      meta.data = meta.data %>% tibble::rownames_to_column(var = colnames(meta.data)[1])
      if (nrow(meta.data) != ncol(count_mat)) {
        message("检测到meta.data的细胞数(nrow)与count_mat的细胞数(ncol)不一致")
        message(paste0("meta.data的细胞数(nrow): ", nrow(meta.data)))
        message(paste0("count_mat的细胞数(ncol): ", ncol(count_mat)))
        stop()
      }
      if (sum(rownames(meta.data) %in% colnames(count_mat)) != ncol(count_mat)) {
        message("检测到meta.data的rownames与count_mat的colnames没有完全一一对应")
        message(paste0("没有在count_mat的colnames中出现的meta.data的rownames: ", paste(head(rownames(meta.data)[!rownames(meta.data) %in% colnames(count_mat)]), collapse = ", ")))
        stop()
      }
      # 自动对齐meta.data和count_mat的细胞顺序
      meta.data = meta.data[colnames(count_mat), ]
    }
  } else {
    meta.data = NULL
  }
  
  ### 合并重复基因名的表达量 ###
  if (gene_id_input) {
    print("-------------- merge duplicated gene-ids --------------")
    if (species == "mouse") {
      gene.table = gene.table.mouse
      gene.table.duplicate = gene.table.mouse.duplicate
      gene.table.remove = gene.table.mouse.remove
    }
    if (species == "human") {
      gene.table = gene.table.human
      gene.table.duplicate = gene.table.human.duplicate
      gene.table.remove = gene.table.human.remove
    }
    gene.table.duplicate$hit = gene.table.duplicate$`gene id` %in% rownames(count_mat)
    gene.table.duplicate$expr = rowSums(count_mat)[gene.table.duplicate$`gene id`]
    gene.table.duplicate$expr[is.na(gene.table.duplicate$expr)] = 0
    library(progress)
    n = length(unique(gene.table.duplicate$`gene name`))
    pb = progress_bar$new(format = "[:bar] :percent :elapsed", total = n)
    for (gene in unique(gene.table.duplicate$`gene name`)) {
      tmp = gene.table.duplicate[gene.table.duplicate$`gene name` %in% gene, ]
      tmp = tmp %>% dplyr::filter(hit == TRUE)
      id.n = length(tmp$hit)
      hit.n = sum(tmp$hit)
      if (hit.n > 1) {
        # print(tmp)
        # print(count_mat[tmp$`gene id`, ][, 1:10])
        sum.expr = colSums(count_mat[tmp$`gene id`, ])
        count_mat[tmp$`gene id`, ] = 0
        count_mat[tmp$`gene id`[1], ] = sum.expr
        # print("⬇⬇⬇⬇⬇")
        # print(count_mat[tmp$`gene id`, ][, 1:10])
        count_mat = count_mat[!(rownames(count_mat) %in% tmp$`gene id`[2:id.n]), ]
      }
      pb$tick()
    }
    rownames(gene.table) = gene.table$`gene id`
    rownames(count_mat) = gene.table[rownames(count_mat), "gene name"]
    message("将多个gene id对1个gene name的情况做加和处理完成!")
  }
  
  ### 2.去除non_protein_coding基因和低检测率基因 ###
  print("-------------- remove non protein_coding genes --------------")
  if (species == "mouse") {
    gene.table = gene.table.mouse
    gene.table.duplicate = gene.table.mouse.duplicate
    gene.table.remove = gene.table.mouse.remove
    genes_non_protein_coding = gene.table.remove$`gene name`
    genes_remove = c(genes_remove, genes_non_protein_coding)
    genes_keep = rownames(count_mat)[!rownames(count_mat) %in% genes_remove]
    count_mat = count_mat[genes_keep, ]
    message(paste0("去除non_protein_coding基因后的count_mat有", ncol(count_mat), "个细胞, ", nrow(count_mat), "个基因"))
  }
  if (species == "human") {
    gene.table = gene.table.human
    gene.table.duplicate = gene.table.human.duplicate
    gene.table.remove = gene.table.human.remove
    genes_non_protein_coding = gene.table.remove$`gene name`
    genes_remove = c(genes_remove, genes_non_protein_coding)
    genes_keep = rownames(count_mat)[!rownames(count_mat) %in% genes_remove]
    count_mat = count_mat[genes_keep, ]
    message(paste0("去除non_protein_coding基因后的count_mat有", ncol(count_mat), "个细胞, ", nrow(count_mat), "个基因"))
  }
  
  # 自动转为稀疏矩阵以节省内存
  count_mat = as.sparse(count_mat)

  ### 3.创建Seurat对象 ###
  if (run_geneQC) {
    print("-------------- QC for gene --------------")
    min.cells = max(c(floor(ncol(count_mat)*0.001), min.cells))
    message(paste0("对于基因, cutoff为:\nmin.cells > ", min.cells))
    seurat = CreateSeuratObject(counts = count_mat, assay = assay, project = project, min.cells = min.cells, min.features = ifelse(run_cellQC, min.features, 0), meta.data = meta.data)
    message(paste0("成功基于上述cutoff, seurat对象有", ncol(seurat), "个细胞, ", nrow(seurat), "个基因"))
  }
  seurat = CreateSeuratObject(counts = count_mat, assay = assay, project = project, min.cells = 0, min.features = ifelse(run_cellQC, min.features, 0), meta.data = meta.data)
  seurat[["orig.ident"]] = project
  
  # ST合法性
  if (isST) {
    # rename cell names
    if (!file.exists(barcode_df)) {
      stop("确保barcode_df合法!")
    }
    barcode.df = read.csv(barcode_df)
    rownames(barcode.df) = paste0(barcode.df$bc_B, barcode.df$bc_A)
    if (sum(colnames(seurat) %in% rownames(barcode.df)) == 0) {
      stop("colnames(seurat)和rownames(barcode.df)完全不对应, 检查前者是否是bBbA, 或者gene x barcode矩阵倒转了!")
    }
    barcode.df = barcode.df[colnames(seurat), ]
    seurat = RenameCells(seurat, new.names = paste0(barcode.df$iB, "x", barcode.df$iA))
    if (!is.null(image_dir)) {
      if (!dir.exists(image_dir)) {
        stop(paste0("请检查image_dir是:", image_dir))
      } else {
        if (numBC == 96) {
          create_tissue_position(96, paste(image_dir, "tissue_positions_list.csv", sep = "/"))
          json = '{"spot_diameter_fullres": 1, "tissue_hires_scalef": 1, "fiducial_diameter_fullres": 96, "tissue_lowres_scalef": 1}'
          cat(json, file = paste(image_dir, "scalefactors_json.json", sep = "/"), fill = FALSE, labels = NULL, append = FALSE)
        } else if (numBC == 50) {
          create_tissue_position(50, paste(image_dir, "tissue_positions_list.csv", sep = "/"))
          json = '{"spot_diameter_fullres": 1, "tissue_hires_scalef": 1, "fiducial_diameter_fullres": 184.32, "tissue_lowres_scalef": 1}'
          cat(json, file = paste(image_dir, "scalefactors_json.json", sep = "/"), fill = FALSE, labels = NULL, append = FALSE)
        }
        image = Read10X_Image(image_dir, image.name = "tissue_lowres_image.png")
        if (sum(rownames(image@coordinates) %in% colnames(seurat)) != ncol(seurat)) {
          message("检测到image@coordinates的rownames与seurat的colnames没有完全一一对应")
          message(paste0("没有在seurat的colnames中出现的image@coordinates的rownames: ", paste(head(rownames(image@coordinates)[!rownames(image@coordinates) %in% colnames(seurat)]), collapse = ", ")))
          stop()
        }
        # 自动对齐image@coordinates和seurat的细胞顺序
        image = image[colnames(seurat)]
      }
    } else {
      stop("请输入image_dir")
    }
  }
  # 添加ST信息
  if (isST) {
    DefaultAssay(image) = assay
    seurat[[project]] = image[colnames(seurat)]
    seurat$bc_B = barcode.df$bc_B
    seurat$bc_A = barcode.df$bc_A
    seurat$bc_full = paste0(barcode.df$bc_B, barcode.df$bc_A)
  }
  # 添加细胞信息
  seurat$cell_id = colnames(seurat)
  
  ### 4.再次QC ###
  if (run_cellQC) {
    print("-------------- QC for cell --------------")
    nCount.min = min.nCount
    nCount.max = floor(quantile(seurat@meta.data[ , paste0("nCount_", assay)], probs = max.nCount.quantile)[[1]])
    message(paste0("对于细胞, cutoff为:\n", nCount.min, " <= nCount <= ", nCount.max))
    nFeature.min = min.features
    nFeature.max = floor(quantile(seurat@meta.data[ , paste0("nFeature_", assay)], probs = max.nFeature.quantile)[[1]])
    message(paste0(nFeature.min, " <= nFeature <= ", nFeature.max))
    cells_valid = seurat@meta.data %>% 
      dplyr::filter(!!as.name(paste0("nCount_", assay)) >= nCount.min & !!as.name(paste0("nCount_", assay)) <= nCount.max & 
                      !!as.name(paste0("nFeature_", assay)) >= nFeature.min & !!as.name(paste0("nFeature_", assay)) <= nFeature.max) %>% 
      dplyr::pull(cell_id)
    seurat = subset(seurat, cells = cells_valid)
    message(paste0("成功基于上述cutoff, seurat对象有", ncol(seurat), "个细胞, ", nrow(seurat), "个基因"))
  }
  
  ### 5.添加额外指标 ###
  # 以下指标含量: gene set的UMI和 / 细胞的总UMI * 100
  # 所以只要细胞的UMI不变(基因表达谱不变), 则对于细胞的subset不会影响原来的指标含量结果
  if (species == "mouse") {
    seurat[["pct.mt"]] = PercentageFeatureSet(seurat, pattern = "^mt-")
    seurat[["pct.ribo"]] = PercentageFeatureSet(seurat, pattern = "^Rp[sl]")
    seurat[["pct.hb"]] = PercentageFeatureSet(seurat, pattern = "^Hb[ab]-")
  }
  if (species == "human") {
    seurat[["pct.mt"]] = PercentageFeatureSet(seurat, pattern = "^MT-")
    seurat[["pct.ribo"]] = PercentageFeatureSet(seurat, pattern = "^RP[SL]")
    seurat[["pct.hb"]] = PercentageFeatureSet(seurat, pattern = "^HB[AB]")
  }
  
  ### 6.计算doublet ###
  if (run_doublet) {
    print("-------------- doublet --------------")
    methods = c("cxds", "bcds", "hybrid", "scDblFinder", "Scrublet", "DoubletDetection", "DoubletFinder")
    seurat = find_doublet(seurat, methods = methods, step = "all")
  }
  
  seurat = RenameCells(seurat, new.names = paste0(cellname_prefix, colnames(seurat), cellname_suffix))
  message("\n")
  
  seurat = NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  
  return(seurat)
}