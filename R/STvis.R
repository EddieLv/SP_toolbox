library(magick)
library(shiny)
library(shinyalert)
library(shinydashboard)
library(png)
library(tools)
library(DT)

library(Seurat)
library(dplyr)
library(tibble)
library(data.table)
library(stringr)

library(RColorBrewer)
library(ggiraph)
library(ggplot2)
library(viridis)
library(grid)

shiny_st = function(seurat, assay = "SCT", slot = "data", image = NULL, python_env = NULL, script = NULL, tooltip = NULL, crop = TRUE) {
  DefaultAssay(seurat) = assay
  
  move.axis.shiny = function(df, x = NULL, y = NULL, numBarcode = NULL, x.num = 0, y.num = 0) {
    # 定义旋转角度和旋转中心坐标
    n.min = (1 - 0.5) * 1080 / numBarcode
    
    df.old = df
    # 赋值时x y对调
    df[, y] = df.old[, y] + y.num*n.min
    df[, x] = df.old[, x] + x.num*n.min
    
    return(df)
  }
  
  shrink.axis.shiny = function(df, x = NULL, y = NULL, numBarcode = NULL, x.scale.factor = 1, y.scale.factor = 1) {
    
    df.old = df
    # 赋值时x y对调
    df[, y] = df.old[, y]*y.scale.factor
    df[, x] = df.old[, x]*x.scale.factor
    
    return(df)
  }
  
  rotate.axis.shiny = function(df, x = NULL, y = NULL, numBarcode = NULL, angle = 0) {
    rad = angle * pi / 180
    
    # 定义旋转角度和旋转中心坐标
    n.min = (1 - 0.5) * 1080 / numBarcode
    n.max = (numBarcode - 0.5) * 1080 / numBarcode
    n.mid = (n.min + n.max) / 2
    center = c(n.mid, n.mid)
    
    df.old = df
    # 赋值时x y对调
    df[, y] = (df.old[, y] - center[2]) * cos(rad) - (df.old[, x] - center[1]) * sin(rad) + center[1] # y1 = y*cos(β) - x*sin(β)
    df[, x] = (df.old[, x] - center[1]) * cos(rad) + (df.old[, y] - center[2]) * sin(rad) + center[2] # x1 = x*cos(β) + y*sin(β)
    
    return(df)
  }
  
  flip.axis.shiny = function(df, x = NULL, y = NULL, numBarcode = NULL, horizontal = F, vertical = F) {
    # (numBarcode + 1 - 0.5 - iA) * 1080 / numBarcode
    n.min = (1 - 0.5) * 1080 / numBarcode
    n.max = (numBarcode - 0.5) * 1080 / numBarcode
    
    df.old = df
    if (horizontal) {
      if (numBarcode == 50) {
        df[, y] = n.max - df.old[, y] + n.min
      }
      if (numBarcode == 96) {
        df[, y] = n.max - df.old[, y] + n.min
      }
    }
    
    if (vertical) {
      if (numBarcode == 50) {
        df[, x] = n.max - df.old[, x] + n.min
      }
      if (numBarcode == 96) {
        df[, x] = n.max - df.old[, x] + n.min
      }
    }
    
    return(df)
  }
  
  create_tissue_position.shiny = function(numBarcode) {
    
    barcode_index = 1:numBarcode
    tissue.positions = data.frame(iB = rep(barcode_index, each = numBarcode),
                                  iA = rep(barcode_index, times = numBarcode),
                                  tissue = 1)
    
    if (!numBarcode %in% c(50, 96)) {
      stop("Please make sure your barcode num is 50 or 96!")
    }
    
    if (numBarcode == 96) {
      # imagecol = (numBarcode + 1 - 0.0 - iA) * 1080 / numBarcode
      # imagerow = (iB - 0.0) * 1080 / numBarcode)
      tissue.positions = tissue.positions %>% mutate(barcodes = paste0(iB, "x", iA),
                                                     imagecol = (numBarcode + 1 - 0.5 - iA) * 1080 / numBarcode,
                                                     imagerow = (iB - 0.5) * 1080 / numBarcode) %>%
        rename(col = iB, row = iA)
    }
    
    if (numBarcode == 50) {
      tissue.positions = tissue.positions %>% mutate(barcodes = paste0(iB, "x", iA),
                                                     imagecol = (iA - 0.5) * 1080 / numBarcode,
                                                     imagerow = (numBarcode + 1 - 0.5 - iB) * 1080 / numBarcode) %>%
        rename(col = iB, row = iA)
    }
    
    rownames(tissue.positions) = tissue.positions$barcodes
    tissue.positions = tissue.positions %>%
      select(c(tissue, col, row, imagecol, imagerow))
    
    return(tissue.positions)
  }
  
  read_spatial = function (numBarcode, spatial_img) {
    image = readPNG(source = spatial_img)
    scale.factors = list(
      "spot_diameter_fullres" = 1,
      "tissue_hires_scalef" = 1,
      "fiducial_diameter_fullres" = ifelse(numBarcode == 50, 184.32, 96),
      "tissue_lowres_scalef" = 1
    )
    
    tissue.positions = create_tissue_position.shiny(numBarcode)
    
    unnormalized.radius = scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
    spot.radius = unnormalized.radius/max(dim(x = image))
    
    return(new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                                                               fiducial = scale.factors$fiducial_diameter_fullres,
                                                                               hires = scale.factors$tissue_hires_scalef, scale.factors$tissue_lowres_scalef),
               coordinates = tissue.positions, spot.radius = spot.radius))
  }
  
  # Build at least n highly-distinct categorical colors. Combines several curated
  # Seurat DiscretePalette sets (~96 distinct) and, if more categories exist, extends
  # with evenly-spaced HCL hues so any number of labels gets a unique color (>=100 ok).
  discrete.palette.shiny = function(n) {
    base = unique(c(
      Seurat::DiscretePalette(36, palette = "polychrome", shuffle = FALSE),
      Seurat::DiscretePalette(32, palette = "glasbey",    shuffle = FALSE),
      Seurat::DiscretePalette(26, palette = "alphabet",   shuffle = FALSE),
      Seurat::DiscretePalette(26, palette = "alphabet2",  shuffle = FALSE)
    ))
    if (n <= length(base)) return(base[seq_len(n)])
    extra = grDevices::hcl(h = seq(15, 375, length.out = n - length(base) + 1)[-1], c = 100, l = 65)
    c(base, extra)
  }

  crop.limits.shiny = function(coordinates, img, crop = TRUE) {
    if (!isTRUE(crop)) {
      return(list(x = c(0, ncol(img)), y = c(0, nrow(img))))
    }

    padded.range = function(values, image.extent) {
      values = values[is.finite(values)]
      if (length(values) == 0) {
        return(c(0, image.extent))
      }

      limits = range(values)
      # Match SpatialPlot(crop = TRUE): frame the spots rather than the full slide,
      # while retaining a small amount of tissue around the outermost spots.
      padding = max(diff(limits) * 0.05, image.extent * 0.005, 1)
      limits = c(limits[1] - padding, limits[2] + padding)
      c(max(0, limits[1]), min(image.extent, limits[2]))
    }

    list(
      x = padded.range(coordinates$x, ncol(img)),
      y = padded.range(coordinates$y, nrow(img))
    )
  }

  make.feature.plot.shiny = function(ann = NULL, anno.df = NULL, alpha = 0.8, pt.size = 0.1, shape = 22, show.feature = NULL, mode = NULL, tooltip = tooltip) {
    annotation = ann[[1]]
    coordinates = ann[[2]]
    img = ann[[3]]
    coordinate.index = match(anno.df$barcode, rownames(coordinates))
    keep = !is.na(coordinate.index)
    coordinates = coordinates[coordinate.index[keep], , drop = FALSE]
    anno.df = anno.df[keep, , drop = FALSE]
    if (nrow(coordinates) == 0) {
      stop("No cells in the current Seurat object have coordinates in image '", image, "'.")
    }
    plot.limits = crop.limits.shiny(coordinates, img, crop = ann[[4]])
    if (is.factor(anno.df[[show.feature]])) {
      coordinates$feature = droplevels(anno.df[[show.feature]])
    } else {
      coordinates$feature = anno.df[[show.feature]]
    }
    
    if (!is.null(tooltip)) {
      coordinates$tooltip = anno.df[[tooltip]]
    }
    
    if (is.numeric(coordinates$feature)) {
      cols = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100)
      coordinates$feature = round(coordinates$feature, 3)
      p = ggplot(coordinates, aes_string(x = "x", y = "y", data_id = "id_stvis", tooltip = ifelse(is.null(tooltip), "feature", "tooltip"))) + 
        annotation +
        geom_point_interactive(aes(fill = feature, alpha = feature), size = pt.size, shape = shape, stroke = NA, color = "black") +
        scale_fill_gradientn(colors = cols) +
        scale_alpha(range = c(alpha, 1)) +
        ylim(plot.limits$y[2], plot.limits$y[1]) + xlim(plot.limits$x[1], plot.limits$x[2]) +
        theme_void() + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = F, clip = "on") +
        theme(legend.position = "top", plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) +
        guides(alpha = "none") +
        labs(fill = show.feature)
    } else {
      n.cat = length(unique(stats::na.omit(as.character(coordinates$feature))))
      cols = discrete.palette.shiny(max(n.cat, 1L))
            
      # only for AI filtering
      if (mode == "on") {
        shape = 0
      }
      
      p = ggplot(coordinates, aes_string(x = "x", y = "y", data_id = "id_stvis", tooltip = ifelse(is.null(tooltip), "feature", "tooltip"))) + 
        annotation +
        geom_point_interactive(aes(fill = feature), size = pt.size, shape = shape, stroke = 0.1, alpha = alpha, color = "black") +
        ylim(plot.limits$y[2], plot.limits$y[1]) + xlim(plot.limits$x[1], plot.limits$x[2]) +
        theme_void() + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = F, clip = "on") +
        theme(legend.spacing.y = unit(0, "cm"), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) +
        scale_fill_manual_interactive(values = cols) +
        guides(alpha = "none", fill = guide_legend(ncol = 1, byrow = T, override.aes = list(size = 5))) +
        labs(fill = show.feature)
    }
    
    return(p)
  }
  
  # Pick the scale factor (hires/lowres) whose scaled coords best fit the stored image.
  # The stored @image resolution varies by loading pipeline: standard Visium often keeps
  # the hires image, Visium HD often keeps the lowres image. Hardcoding one mis-scales the
  # other (e.g. HD with hires pushes spots ~10x off-canvas, so no spots are visible).
  pick.scale.shiny = function(img_obj, max.x, max.y) {
    img = img_obj@image
    iw = ncol(img); ih = nrow(img)
    sfs = img_obj@scale.factors
    cand = c("hires", "lowres")
    score = sapply(cand, function(s) {
      fill = max(max.x * sfs[[s]] / iw, max.y * sfs[[s]] / ih)
      if (fill > 1.05) (fill - 1) * 1000 else (1 - fill)  # penalise overflow, prefer best fill
    })
    sfs[[ cand[which.min(score)] ]]
  }

  add_image = function (seurat, image.name, crop = TRUE) {
    img_obj = seurat@images[[image.name]]
    img = img_obj@image

    if (inherits(img_obj, "VisiumV2")) {
      # VisiumV2 (Seurat v5): coords are full-resolution pixel space (x = column, y = row).
      # Scale to the stored image and map directly to the plot axes (no rotate/flip): this
      # aligns spots to the image for both square (Visium) and rectangular (HD) tissues.
      coordinates = GetTissueCoordinates(seurat, image = image.name, scale = NULL)[, 1:2]
      cells.use = intersect(colnames(seurat), rownames(coordinates))
      coordinates = coordinates[cells.use, , drop = FALSE]
      colnames(coordinates) = c("x", "y")
      sf = pick.scale.shiny(img_obj, max(coordinates$x), max(coordinates$y))
      coordinates$x = coordinates$x * sf
      coordinates$y = coordinates$y * sf
      coordinates$id_stvis = seurat@meta.data[rownames(coordinates), "id_stvis"]
    } else {
      # VisiumV1 (legacy DBIT-seq / Seurat v4): raw grid coords scaled by lowres, then
      # transposed (rotate 90 + flip horizontal) to match the image orientation.
      if (utils::packageVersion("Seurat") >= "5.0.0") {
        coordinates = GetTissueCoordinates(seurat, image = image.name, scale = NULL)[, 1:2]
      } else {
        coordinates = GetTissueCoordinates(seurat, image = image.name)[, 1:2]
      }
      cells.use = intersect(colnames(seurat), rownames(coordinates))
      coordinates = coordinates[cells.use, , drop = FALSE]
      colnames(coordinates) = c("x", "y")
      coordinates = coordinates %>%
        mutate(x = x * img_obj@scale.factors$lowres,
               y = y * img_obj@scale.factors$lowres)
      coordinates = rotate.axis.shiny(coordinates, x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), angle = 90)
      coordinates = flip.axis.shiny(coordinates, x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), horizontal = T)
      coordinates$id_stvis = seurat@meta.data[rownames(coordinates), "id_stvis"]
    }

    img_grob = grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
    # ggplot2 v4.0+ 严格按数据范围裁剪 annotation_custom，ymax 必须为正值才在 ylim(nrow,0) 范围内
    annotation = annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = nrow(img))

    return(list(annotation, coordinates, img, crop))
  }
  
  ai.filter = function(image = NULL, df = NULL, python_env = NULL, script = NULL, barcodeNum = NULL, thre = NULL, prefix = NULL) {
    writePNG(image, target = paste(getwd(), "STvis.AI.filter.genger.temp.png", sep = "/"))
    system(paste(python_env, script,
                 paste0("--img=", paste(getwd(), "STvis.AI.filter.genger.temp.png", sep = "/")),
                 paste0("--barcodeNum=", barcodeNum),
                 paste0("--thre=", thre), sep = " "))
    ai.filtered.pixels = read.csv(paste("/home/biogenger/Downloads/STvis.AI.filter.genger.filterd.pixels.csv", sep = "/"), header = F)
    system(paste("rm -rf", paste(getwd(), "STvis.AI.filter.genger.temp.png", sep = "/")))
    system(paste("rm -rf", paste(getwd(), "/home/biogenger/Downloads/STvis.AI.filter.genger.filterd.pixels.csv", sep = "/")))
    colnames(ai.filtered.pixels) = "barcode"
    # add important feature!
    df[["ai.filter"]][df[["barcode"]] %in% paste0(prefix, ai.filtered.pixels$barcode)] = "filtered"
    
    return(df[["ai.filter"]])
  }
  
  
  ui = dashboardPage(dashboardHeader(title = tags$div(style = "white-space: pre-wrap; word-wrap: break-word; line-height: 1.2;",
                                                      "ST Visualization Tool (STvis v1.0.0)\n---please watch the [Instructions] first---"), titleWidth = 450),
                     dashboardSidebar(width = 600,
                                      actionButton(inputId = "info", label = "Instructions"), shiny::hr(),
                                      
                                      fluidRow(
                                        tags$div(
                                          style = "display: flex; align-items: center;",
                                          sliderInput(inputId = "thre", label = "threshold for filtering [0-100]", min = 0, max = 100, value = 50, step = 1, width = "60%"),
                                          actionButton(inputId = "ai.ok", label = "Power BY AI", width = "20%"),
                                          selectInput(inputId = "filter.mode", label = "Filter mode", choices = c("on", "off"), selected = "off", width = "20%"))), shiny::hr(),
                                      
                                      fluidRow(
                                        column(width = 4,
                                               selectInput(inputId = "sampleInput", label = "Select sample", choices = NULL, selected = NULL, width = "100%")),
                                        column(width = 4,
                                               selectInput(inputId = "shapeInput", label = "Select shape", choices = c(22), selected = 22, width = "100%")),
                                        column(width = 4,
                                               selectInput(inputId = "featureInput", label = "Select feature", choices = c("orig.ident"), selected = "orig.ident", width = "100%"))), shiny::hr(),
                                      
                                      fluidRow(
                                        tags$div(
                                          style = "display: flex; align-items: center;",
                                          selectInput(inputId = "subsetFeature", label = "Select feature to subset", choices = c("orig.ident"), selected = "orig.ident", width = "40%"),
                                          textInput(inputId = "subsetString", label = "Select specific features based on selected feature.\ne.g. [typeA,typeB]", value = "", width = "40%"),
                                          actionButton(inputId = "subset.ok", label = "Select", width = "10%")),
                                        tags$div(
                                          style = "display: flex; align-items: center;",
                                          actionButton(inputId = "recover.ok", label = "Back to all celltypes!", width = "50%"))), shiny::hr(),
                                      
                                      fluidRow(
                                        column(width = 6,
                                               sliderInput(inputId = "alphaValue", label = "Spot.alpha [0-1]", min = 0, max = 1, value = 0.8, step = 0.01)),
                                        column(width = 6,
                                               sliderInput(inputId = "spotSize", label = "Spot.size [0-1]", min = 0, max = 1, value = 0.1, step = 0.01))), shiny::hr(),
                                      
                                      fluidRow(
                                        column(width = 6,
                                               actionButton(inputId = "vertical", label = "Flip by vertical", width = "85%"),
                                               actionButton(inputId = "horizontal", label = "Flip by horizontal", width = "85%")),
                                        tags$div(
                                          style = "display: flex; align-items: center;",
                                          numericInput(inputId = "angle", label = "Rotate by angle [-360, 360, 1]", value = 0, min = -360, max = 360, step = 1, width = "70%"),
                                          actionButton(inputId = "angle.ok", label = "√", width = "10%"))),
                                      
                                      tags$div(
                                        tags$div(
                                          style = "flex-basis: 50%; margin-bottom: 10px;",
                                          tags$div(
                                            style = "display: flex; justify-content: space-between;",
                                            tags$div(
                                              style = "display: flex; align-items: center;",
                                              numericInput(inputId = "xMove", label = "Move spots horizontally [-96, 96, 1]", value = 0, min = -96, max = 96, step = 1, width = "90%"),
                                              actionButton(inputId = "xMove.ok", label = "√", width = "10%")
                                            ),
                                            tags$div(
                                              style = "display: flex; align-items: center;",
                                              numericInput(inputId = "yMove", label = "Move spots vertically [-96, 96, 1]", value = 0, min = -96, max = 96, step = 1, width = "90%"),
                                              actionButton(inputId = "yMove.ok", label = "√", width = "10%")
                                            )
                                          )),
                                        
                                        tags$div(
                                          style = "flex-basis: 50%;",
                                          tags$div(
                                            style = "display: flex; justify-content: space-between;",
                                            tags$div(
                                              style = "display: flex; align-items: center;",
                                              numericInput(inputId = "xShrink", label = "Shrink spots horizontally [0, 5, 0.1]", value = 0, min = -0, max = 5, step = 0.1, width = "90%"),
                                              actionButton(inputId = "xShrink.ok", label = "√", width = "10%")
                                            ),
                                            tags$div(
                                              style = "display: flex; align-items: center;",
                                              numericInput(inputId = "yShrink", label = "Shrink spots vertically [0, 5, 0.1]", value = 0, min = -0, max = 5, step = 0.1, width = "90%"),
                                              actionButton(inputId = "yShrink.ok", label = "√", width = "10%")
                                            )
                                          ))
                                      ), shiny::hr(),
                                      
                                      tags$div(
                                        style = "display: flex; align-items: center;",
                                        textInput(inputId = "geneInput", label = "Select one gene to visualize", value = ""),
                                        actionButton(inputId = "gene.ok", label = "√")
                                      ), shiny::hr(),
                                      
                                      tags$div(
                                        style = "display: flex; align-items: center;",
                                        textInput(inputId = "labelInput", label = "Set label for selected spots", value = "", width = "50%"),
                                        actionButton(inputId = "confirm", label = "Confirm")), shiny::hr(),
                                      
                                      tags$div(
                                        column(width = 6,
                                               actionButton(inputId = "stopApp", label = "Quit")))),
                     dashboardBody(h4(strong("===Attention==="), br(),
                                      "Do not change any attributes before you confirmed selected spots!", br(),
                                      "Otherwise you will lose your selected spots!", style = "color:red;"),
                                   girafeOutput("Plot1", width = "100%", height = "1280px")),
                     tags$head(
                       tags$style(HTML(".sidebar {width: 600px;}
                                 .my-row .form-group {margin-top: 0; margin-bottom: 0;}
                                 .content-wrapper {margin-left: 0;}
                                 .input-group {margin-right: 0;}")))
  )
  
  server = function(input, output, session) {
    options(shiny.maxRequestSize = 100*1024^2)
    
    if (is.null(image)) {
      stop("Please set image!")
    }

    if (!image %in% Images(seurat)) {
      stop("Image '", image, "' was not found. Available images: ", paste(Images(seurat), collapse = ", "))
    }

    image.coordinates = GetTissueCoordinates(seurat, image = image)
    image.cells = intersect(colnames(seurat), rownames(image.coordinates))
    if (length(image.cells) == 0) {
      stop("No cells in the current Seurat object have coordinates in image '", image, "'.")
    }

    if (length(Images(seurat)) > 1 || length(image.cells) < ncol(seurat)) {
      if (length(Images(seurat)) > 1) {
        message(paste("Detect", length(Images(seurat)), "images; using", image))
      }
      seurat.backup = seurat
      seurat = subset(seurat.backup, cells = image.cells, slot = slot)
      seurat@meta.data = seurat@meta.data[ , colnames(seurat.backup@meta.data)]
      img = list(seurat@images[[image]]); names(img) = image
      seurat@images = img
    }
    
    # if (sum(str_detect(colnames(seurat), "x")) < dim(seurat)[2]) {
    #   stop("Your cell id_stvis must be formatted like 1x1, 1x2 or sample_1sx1, sample_1x2 ...")
    # }
    
    sampleChoice = as.character(unique(seurat$orig.ident))
    shapeChoice = c(22, 21)
    featureChoice = colnames(seurat@meta.data)
    # for (i in featureChoice) {
    #   seurat@meta.data[is.na(seurat[[i]]), i] = "STvis"
    #   seurat@meta.data[is.null(seurat[[i]]), i] = "STvis"
    # }
    
    # clean.barcodes = str_match(colnames(seurat), pattern = "\\d+x\\d+")[ , 1]
    image.coordinates = GetTissueCoordinates(seurat, image = image)
    image.coordinates = image.coordinates[rownames(seurat@meta.data), , drop = FALSE]
    clean.barcodes = paste0(image.coordinates[, 1], "x", image.coordinates[, 2])
    # prefix = gsub(clean.barcodes[1], "", colnames(seurat)[1])
    prefix = seurat$orig.ident[1]
    # add important feature!
    seurat$barcodeB_stvis = str_split(clean.barcodes, "x", simplify = T)[ , 1]
    seurat$barcodeA_stvis = str_split(clean.barcodes, "x", simplify = T)[ , 2]
    seurat$id_stvis = 1:dim(seurat)[2]
    
    updateSelectInput(session, inputId = "sampleInput", label = "Select sample", choices = sampleChoice, selected = sampleChoice[1])
    updateSelectInput(session, inputId = "shapeInput", label = "Select shape", choices = shapeChoice, selected = 22)
    updateSelectInput(session, inputId = "featureInput", label = "Select feature", choices = featureChoice, selected = "orig.ident")
    updateSelectInput(session, inputId = "subsetFeature", label = "Select feature to subset", choices = featureChoice, selected = "orig.ident")
    
    df = reactiveValues()
    df.backup = reactiveValues()
    for (i in c(featureChoice, "id_stvis")) {
      if (is.factor(seurat@meta.data[[i]])) {
        df[[i]] = droplevels(seurat@meta.data[[i]])
        df.backup[[i]] = droplevels(seurat@meta.data[[i]])
      } else {
        df[[i]] = seurat@meta.data[[i]]
        df.backup[[i]] = seurat@meta.data[[i]]
      }
      df$barcode = rownames(seurat@meta.data)
      df.backup$barcode = rownames(seurat@meta.data)
    }
    
    rv = reactiveValues(sNr = "1", ann = NULL)
    observeEvent(input$sampleInput, {
      rv$sNr = sampleChoice[1]
      rv$ann = add_image(seurat, image.name = image, crop = crop)
    })
    
    observeEvent(input$gene.ok, {
      mat.df = GetAssayData(seurat, assay = assay, slot = slot)
      if (!input$geneInput %in% rownames(mat.df)) {
        shinyalert(paste(input$geneInput, "does not exists in", Images(seurat)[1], "-", assay, "-", slot, "!"))
        return(NULL)
      }
      mat.df = FetchData(seurat, vars = input$geneInput, slot = slot) %>% rownames_to_column(var = "barcode")
      df.backup[[input$geneInput]] = mat.df[[input$geneInput]][mat.df$barcode %in% df.backup$barcode]
      df[[input$geneInput]] = mat.df[[input$geneInput]][mat.df$barcode %in% df$barcode]
      updateSelectInput(session, inputId = "featureInput", label = "Select feature", choices = names(df), selected = input$geneInput)
      updateSelectInput(session, inputId = "subsetFeature", label = "Select feature to subset", choices = names(df), selected = "orig.ident")
    })
    
    cols = reactiveValues(selected.fill = "transparent")
    observeEvent(input$subset.ok, {
      subset.features = strsplit(input$subsetString, split = ",")[[1]]
      
      if (all(subset.features %in% df[[input$subsetFeature]])) {
        for (i in names(df)[names(df) != input$subsetFeature]) {
          df[[i]] = df[[i]][df[[input$subsetFeature]] %in% subset.features]
        }
        
        if (is.factor(df[[input$subsetFeature]])) {
          df[[input$subsetFeature]] = droplevels(df[[input$subsetFeature]])
        }
        df[[input$subsetFeature]] = df[[input$subsetFeature]][df[[input$subsetFeature]] %in% subset.features]
        
      } else {
        
        for (i in names(df)[names(df) != input$subsetFeature]) {
          df[[i]] = df.backup[[i]][df.backup[[input$subsetFeature]] %in% subset.features]
        }
        df[[input$subsetFeature]] = df.backup[[input$subsetFeature]][df.backup[[input$subsetFeature]] %in% subset.features]
        
      }
      
    })
    
    observeEvent(input$recover.ok, {
      for (i in names(df)) {
        df[[i]] = df.backup[[i]]
      }
    })
    
    observeEvent(input$vertical, {
      rv$ann[[2]] = flip.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), horizontal = T)
    })
    
    observeEvent(input$horizontal, {
      rv$ann[[2]] = flip.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), vertical = T)
    })
    
    observeEvent(input$angle.ok, {
      rv$ann[[2]] = rotate.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), angle = input$angle)
    })
    
    observeEvent(input$xMove.ok, {
      rv$ann[[2]] = move.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), x.num = input$xMove)
    })
    
    observeEvent(input$yMove.ok, {
      rv$ann[[2]] = move.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), y.num = input$yMove)
    })
    
    observeEvent(input$xShrink.ok, {
      rv$ann[[2]] = shrink.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), x.scale.factor = input$xShrink)
    })
    
    observeEvent(input$yShrink.ok, {
      rv$ann[[2]] = shrink.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), y.scale.factor = input$yShrink)
    })
    
    output$Plot1 = renderGirafe({
      withProgress(message = "Updating plot", value = 0,
                   {
                     df.tmp = as.data.frame(reactiveValuesToList(df))
                     x = girafe(ggobj = make.feature.plot.shiny(ann = rv$ann, anno.df = df.tmp, alpha = input$alphaValue, pt.size = input$spotSize*10, mode = input$filter.mode,
                                                                shape = as.integer(input$shapeInput), show.feature = ifelse(is.null(input$featureInput), "orig.ident", input$featureInput), 
                                                                tooltip = tooltip),
                                width_svg = 12, height_svg = 10)
                     x = girafe_options(x,
                                        opts_zoom(min = 1, max = 10),
                                        opts_tooltip(opacity = 1),
                                        opts_toolbar(position = "topright",
                                                     saveaspng = TRUE,
                                                     pngname = Images(seurat)[1],
                                                     list(lasso_select = "lasso",
                                                          zoom_on = "zoom",
                                                          zoom_reset = "recover",
                                                          saveaspng = "download png"),
                                                     hidden = c("lasso_deselect", "zoom_rect")),
                                        opts_selection(type = "multiple", css = paste0("fill:", cols$selected.fill, ";stroke:transparent;stroke-width:0.1px")),
                                        opts_selection_key(css = "stroke:black;r:1pt;"),
                                        opts_hover(css = "fill:wheat;stroke:black;stroke-width:1px;cursor:pointer;"),
                                        opts_hover_key(css = "stroke:black;r:5pt;cursor:pointer;"))
                     x
                   })
    })
    
    observeEvent(input$ai.ok, {
      df[["ai.filter"]] = rep("exist", length(df[["id_stvis"]]))
      df.backup[["ai.filter"]] = rep("exist", length(df.backup[["id_stvis"]]))
      featureChoice = unique(c(featureChoice, "ai.filter"))
      updateSelectInput(session, inputId = "featureInput", label = "Select feature", choices = featureChoice, selected = "ai.filter")
      updateSelectInput(session, inputId = "subsetFeature", label = "Select feature to subset", choices = featureChoice, selected = "ai.filter")
      
      df[["ai.filter"]] = ai.filter(rv$ann[[3]], df, python_env = python_env, script = script, barcodeNum = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), thre = input$thre, prefix = prefix)
      df.backup[["ai.filter"]] = ai.filter(rv$ann[[3]], df.backup, python_env = python_env, script = script, barcodeNum = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), thre = input$thre, prefix = prefix)
      
      session$sendCustomMessage(type = "Plot1_set", message = character(0))
    })
    
    observeEvent(input$confirm, {
      shinyalert("Remember to click [quit] to save your confirmed labels!")
      ids.selected = as.numeric(input$Plot1_selected)
      if (is.factor(df[[input$featureInput]])) {
        levs = levels(df[[input$featureInput]])
        df[[input$featureInput]] = as.character(df[[input$featureInput]])
        df[[input$featureInput]][which(df$id_stvis %in% ids.selected)] = input$labelInput
        df[[input$featureInput]] = factor(df[[input$featureInput]], levels = unique(c(levs, input$labelInput)))
        
        levs = levels(df.backup[[input$featureInput]])
        df.backup[[input$featureInput]] = as.character(df.backup[[input$featureInput]])
        df.backup[[input$featureInput]][which(df.backup$id_stvis %in% ids.selected)] = input$labelInput
        df.backup[[input$featureInput]] = factor(df.backup[[input$featureInput]], levels = unique(c(levs, input$labelInput)))
      } else {
        df[[input$featureInput]][which(df$id_stvis %in% ids.selected)] = input$labelInput
        
        df.backup[[input$featureInput]][which(df.backup$id_stvis %in% ids.selected)] = input$labelInput
      }
      session$sendCustomMessage(type = "Plot1_set", message = character(0))
    })
    
    observe({
      
      if (input$stopApp > 0) {
        spots.df = df.backup$barcode
        spots.seurat = rownames(seurat@meta.data)
        
        featureChoice = featureChoice[!featureChoice %in% c("barcodeB_stvis", "barcodeA_stvis", "id_stvis")]
        # this function does not change the raw columns of seurat but add another renamed column!
        for (i in featureChoice) {
          if (is.factor(seurat@meta.data[[i]])) {
            levs = levels(droplevels(seurat@meta.data[[i]]))
            seurat@meta.data[[i]] = as.character(seurat@meta.data[[i]])
            seurat@meta.data[[i]][spots.seurat %in% spots.df] = as.character(df.backup[[i]])
            seurat@meta.data[[i]] = factor(seurat@meta.data[[i]], levels = c(levs, levels(df.backup[[i]])[!levels(df.backup[[i]]) %in% levs]))
          } else {
            seurat@meta.data[[i]][spots.seurat %in% spots.df] = df.backup[[i]]
          }
        }

        if ("ai.filter" %in% names(df.backup)) {
          seurat@meta.data[["ai.filter"]][spots.seurat %in% spots.df] = df.backup[["ai.filter"]]
        }
        
        if (exists("seurat.backup")) {
          coordinates = rv$ann[[2]][ , c("x", "y")] %>% 
            rotate.axis.shiny(x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), angle = 90) %>% 
            flip.axis.shiny(x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), horizontal = T)

          cols.backup = colnames(seurat.backup@meta.data)[!colnames(seurat.backup@meta.data) %in% c("barcodeB_stvis", "barcodeA_stvis", "id_stvis")]
          for (i in cols.backup) {
            if (is.factor(seurat@meta.data[[i]])) {
              levs = levels(droplevels(seurat.backup@meta.data[[i]]))
              seurat.backup@meta.data[[i]] = as.character(seurat.backup@meta.data[[i]])
              seurat.backup@meta.data[[i]][rownames(seurat.backup@meta.data) %in% rownames(seurat@meta.data)] = as.character(seurat@meta.data[[i]])
              seurat.backup@meta.data[[i]] = factor(seurat.backup@meta.data[[i]], levels = c(levs, levels(seurat@meta.data[ , i])[!levels(seurat@meta.data[ , i]) %in% levs]))
            } else {
              seurat.backup@meta.data[rownames(seurat@meta.data), i] = seurat@meta.data[ , i]
            }
          }
          
          if ("ai.filter" %in% names(df.backup)) {
            seurat.backup@meta.data[["ai.filter"]] = "exist"
            seurat.backup@meta.data[rownames(seurat@meta.data), "ai.filter"] = seurat@meta.data[ , "ai.filter"]
          }
          
          seurat.backup@meta.data = seurat.backup@meta.data[ , !(colnames(seurat.backup@meta.data) %in% c("barcodeB_stvis", "barcodeA_stvis", "id_stvis"))]
          stopApp(returnValue = seurat.backup)
          
        } else {
          coordinates = rv$ann[[2]][ , c("x", "y")] %>% 
            rotate.axis.shiny(x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), angle = 90) %>% 
            flip.axis.shiny(x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB_stvis) > 50, 96, 50), horizontal = T)
          
          seurat@meta.data = seurat@meta.data[ , !(colnames(seurat@meta.data) %in% c("barcodeB_stvis", "barcodeA_stvis", "id_stvis"))]
          stopApp(returnValue = seurat)
          
        }
        
      }
      
    })
    
    observeEvent(input$info, {
      showModal(modalDialog(title = "Instructions",
                            HTML("* Copyright (c) 2023.05.12, genger<br>* All rights reserved.<br>",
                                 "----------------------------------------------------------------------------------------------------------------<br>",
                                 "If encountered with bug, please contact with wechat <13958598285><br>",
                                 "Attention: currently only supports for one sample!<br>",
                                 "----------------------------------------------------------------------------------------------------------------<br>",
                                 "*帮助文档请查看https://github.com/EddieLv/STvis, 记得star一下我:)"),
                            size = "l",
                            easyClose = TRUE,
                            footer = NULL))
    })
    
  }
  
  runApp(list(ui = ui, server = server), launch.browser = getOption("/usr/bin/google-chrome", interactive()))
}

