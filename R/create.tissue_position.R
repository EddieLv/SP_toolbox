library(dplyr)

### create tissue_positions_list.csv with specific number of barcodes ###
create_tissue_position = function(numBarcode, out_path) {
  
  barcode_index = 1:numBarcode
  tissue.positions = data.frame(iB = rep(barcode_index, each = numBarcode), iA = rep(barcode_index, times = numBarcode), tissue = 1)
  
  if (!numBarcode %in% c(50, 96)) {
    stop("Please make sure your barcode num is 50 or 96!")
  }
  
  if (numBarcode == 96) {
    # imagecol = (numBarcode + 1 - 0.0 - iA) * 1080 / numBarcode
    # imagerow = (iB - 0.0) * 1080 / numBarcode)
    tissue.positions = tissue.positions %>% dplyr::mutate(barcodes = paste0(iB, "x", iA), imagecol = (numBarcode + 1 - 0.5 - iA) * 1080 / numBarcode, imagerow = (iB - 0.5) * 1080 / numBarcode) %>%
      dplyr::rename(col = iB, row = iA)
  }
  
  if (numBarcode == 50) {
    tissue.positions = tissue.positions %>% dplyr::mutate(barcodes = paste0(iB, "x", iA), imagecol = (iA - 0.5) * 1080 / numBarcode, imagerow = (numBarcode + 1 - 0.5 - iB) * 1080 / numBarcode) %>% 
      dplyr::rename(col = iB, row = iA)
  }
  
  rownames(tissue.positions) = tissue.positions$barcodes
  tissue.positions = tissue.positions %>% dplyr::select(c(tissue, col, row, imagecol, imagerow))
  write.table(tissue.positions, file = out_path, col.names = F, sep = ",")
  
}
