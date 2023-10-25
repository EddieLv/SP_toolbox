library(magick)
library(shiny)
library(shinyalert)
library(shinydashboard)
library(png)
library(tools)
library(DT)
library(ggnewscale)

library(Seurat)
library(dplyr)
library(tibble)
library(data.table)
library(stringr)

library(ggiraph)
library(ggplot2)
library(viridis)
library(grid)

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

create.trian.coord = function(pairs, type = "up", shrink = 1) {
  if (type == "up") {
    x = c(0, 0, 1) * shrink
    y = c(0, 1, 1) * shrink
  }
  if (type == "low") {
    x = c(0, 1, 1) * shrink
    y = c(0, 0, 1) * shrink
  }

  mat = do.call(
    rbind,
    apply(pairs, 1, function (row) {
      a = row[1]
      b = row[2]
      data.frame(
        x = x + a,
        y = y + b,
        group = paste(a, b, sep = "-")
      )
    }))

  return(mat)
}

Spatial2Featureplot_genger = function(seurat, orig.ident = NULL, image = NULL, show.image = F, shrink = 1, featureA = NULL, featureB = NULL, colA = NULL, colB = NULL, show.label = T, theme.dark = F) {
  meta.df = seurat@meta.data[seurat$orig.ident == orig.ident, ]

  numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50)

  if (is.null(colA)) {
    if (theme.dark) {
      colA = viridis(100, option = "D")
    } else {
      colA = RColorBrewer::brewer.pal(9, "Reds")
    }
  }
  if (is.null(colB)) {
    if (theme.dark) {
      colB = viridis(100, option = "F")
    } else {
      colB = RColorBrewer::brewer.pal(9, "Blues")
    }
  }

  if (theme.dark) {
    boarder.col = "white"
  } else {
    boarder.col = "grey"
  }

  pairs = merge(1:numBarcode, 1:numBarcode)
  # get up-triangle coordinates
  upper = create.trian.coord(pairs, type = "up", shrink = shrink)
  colnames(upper) = c(paste0("up.", colnames(upper)[1:2]), "group")
  # get down-triangle coordinates
  lower = create.trian.coord(pairs, type = "low", shrink = shrink)[1:2]
  colnames(lower) = paste0("low.", colnames(lower))
  # combine
  upper_lower = cbind(upper, lower)
  upper_lower$cell = gsub("-", "x", upper_lower$group)
  upper_lower$cell.up = paste0(upper_lower$up.x, "x", upper_lower$up.y)
  upper_lower$cell.low = paste0(upper_lower$low.x, "x", upper_lower$low.y)

  # add spatial coordinates
  coordinates = data.frame(iB = rep(seq(1, numBarcode + 1, 1*shrink), each = length(seq(1, numBarcode + 1, 1*0.1))),
                           iA = rep(seq(1, numBarcode + 1, 1*shrink), times = length(seq(1, numBarcode + 1, 1*0.1)))) %>%
    mutate(imagerow = (numBarcode + 1 - iB) * 1080 / numBarcode,
           imagecol = (iA - 1) * 1080 / numBarcode)
  coordinates$cell = paste0(coordinates$iB, "x", coordinates$iA)
  coordinates[ , c("imagerow", "imagecol")] = coordinates[ , c("imagerow", "imagecol")] %>%
    mutate(imagerow = imagerow * seurat@images[[image]]@scale.factors$lowres,
           imagecol = imagecol * seurat@images[[image]]@scale.factors$lowres) %>%
    rotate.axis.shiny(x = "imagerow", y = "imagecol", numBarcode = numBarcode, angle = 90) %>%
    flip.axis.shiny(x = "imagerow", y = "imagecol", numBarcode = numBarcode, horizontal = T)
  imagre.row = coordinates$imagerow; names(imagre.row) = coordinates$cell
  imagre.col = coordinates$imagecol; names(imagre.col) = coordinates$cell

  upper_lower$up.x.image = unname(imagre.row[upper_lower$cell.up])
  upper_lower$up.y.image = unname(imagre.col[upper_lower$cell.up])
  upper_lower$low.x.image = unname(imagre.row[upper_lower$cell.low])
  upper_lower$low.y.image = unname(imagre.col[upper_lower$cell.low])

  clean.barcodes = str_match(rownames(meta.df), pattern = "\\d+x\\d+")[ , 1]
  prefix = gsub(clean.barcodes[1], "", rownames(meta.df)[1])
  df = meta.df %>%
    rownames_to_column(var = "barcode") %>%
    mutate("cell" = gsub(prefix, "", barcode))

  df.res = merge(upper_lower, df, by = "cell", all.x = T) %>% na.omit()

  img = seurat@images[[image]]@image
  img_grob = grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))

  if (show.image) {
    p = ggplot(df.res) +
      annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
      geom_polygon(aes(up.x.image, up.y.image, fill = !!as.name(featureA), group = group), colour = boarder.col, linewidth = shrink*0.5) +
      scale_fill_gradientn(colors = colA) +
      new_scale_fill() +
      geom_polygon(aes(low.x.image, low.y.image, fill = !!as.name(featureB), group = group), colour = boarder.col, linewidth = shrink*0.5) +
      scale_fill_gradientn(colours = colB) +
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = T, clip = "on") +
      ylim(nrow(img), 0) + xlim(0, ncol(img)) +
      theme_void() +
      theme(legend.position = "top",
            aspect.ratio = 1)
  } else {
    df.res[ , c("up.x.image", "up.y.image")] = df.res[ , c("up.x.image", "up.y.image")] %>% flip.axis.shiny(x = "up.x.image", y = "up.y.image", numBarcode = numBarcode, horizontal = T)
    df.res[ , c("low.x.image", "low.y.image")] = df.res[ , c("low.x.image", "low.y.image")] %>% flip.axis.shiny(x = "low.x.image", y = "low.y.image", numBarcode = numBarcode, horizontal = T)
    p = ggplot(df.res) +
      geom_polygon(aes(up.x.image, up.y.image, fill = !!as.name(featureA), group = group), colour = boarder.col, linewidth = shrink*0.5) +
      scale_fill_gradientn(colors = colA) +
      new_scale_fill() +
      geom_polygon(aes(low.x.image, low.y.image, fill = !!as.name(featureB), group = group), colour = boarder.col, linewidth = shrink*0.5) +
      scale_fill_gradientn(colours = colB) +
      scale_x_continuous(expand = c(0, 0), breaks = sort(unique(coordinates$imagecol)) + 5.625, labels = 1:(numBarcode + 1), limits = c(0, 1080), sec.axis = dup_axis()) +
      scale_y_continuous(expand = c(0, 0), breaks = sort(unique(coordinates$imagecol)) + 5.625, labels = 1:(numBarcode + 1), limits = c(0, 1080), sec.axis = dup_axis()) +
      theme_void() +
      theme(legend.position = "top",
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
            axis.text.y.left = element_blank(),
            axis.ticks.y.left = element_blank())
  }

  if (theme.dark) {
    p = p + theme(panel.background = element_rect(fill = "black"))
  }

  return(p)
}

shiny_st = function(seurat, assay = "SCT", slot = "data", image = NULL, python_env = NULL, script = NULL) {
  DefaultAssay(seurat) = assay

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

  make.feature.plot.shiny = function(ann = NULL, anno.df = NULL, alpha = 0.8, pt.size = 0.1, shape = 22, show.feature = NULL, mode = NULL) {
    annotation = ann[[1]]
    coordinates = ann[[2]]
    img = ann[[3]]
    coordinates = coordinates[anno.df$barcode, ]
    if (is.factor(anno.df[[show.feature]])) {
      coordinates$feature = droplevels(anno.df[[show.feature]])
    } else {
      coordinates$feature = anno.df[[show.feature]]
    }

    if (is.numeric(coordinates$feature)) {
      cols = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100)
      p = ggplot(coordinates, aes(x = x, y = y, data_id = id, tooltip = round(feature, 3))) +
        annotation +
        geom_point_interactive(aes(fill = feature, alpha = feature), size = pt.size, shape = shape, stroke = NA, color = "black") +
        scale_fill_gradientn(colors = cols) +
        scale_alpha(range = c(alpha, 1)) +
        ylim(nrow(img), 0) + xlim(0, ncol(img)) +
        theme_void() + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = F, clip = "on") +
        theme(aspect.ratio = 1, legend.position = "top", plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) +
        guides(alpha = "none") +
        labs(fill = show.feature)
    } else {
      cols = c("#F6222EFF", "#5A5156FF", "#E4E1E3FF", "#FE00FAFF", "#16FF32FF", "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", "#2ED9FFFF", "#DEA0FDFF", "#AA0DFEFF",
               "#F8A19FFF", "#325A9BFF", "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", "#782AB6FF",
               "#AAF400FF", "#BDCDFFFF", "#822E1CFF", "#B5EFB5FF", "#7ED7D1FF", "#1C7F93FF", "#D85FF7FF", "#683B79FF", "#66B0FFFF", "#3B00FBFF")
      # only for AI filtering
      if (mode == "on") {
        shape = 0
      }

      p = ggplot(coordinates, aes(x = x, y = y, data_id = id, tooltip = feature)) +
        annotation +
        geom_point_interactive(aes(fill = feature), size = pt.size, shape = shape, stroke = 0.1, alpha = alpha, color = "black") +
        ylim(nrow(img), 0) + xlim(0, ncol(img)) +
        theme_void() + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = F, clip = "on") +
        theme(aspect.ratio = 1, legend.spacing.y = unit(0, "cm"), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) +
        scale_fill_manual_interactive(values = cols) +
        guides(alpha = "none", fill = guide_legend(ncol = 1, byrow = T, override.aes = list(size = 5))) +
        labs(fill = show.feature)
    }

    return(p)
  }

  add_image = function (seurat) {

    image = Images(seurat)[1]
    coordinates = GetTissueCoordinates(seurat, image = image) %>%
      mutate(x = imagerow * seurat@images[[image]]@scale.factors$lowres,
             y = imagecol * seurat@images[[image]]@scale.factors$lowres)
    coordinates = rotate.axis.shiny(coordinates, x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), angle = 90)
    coordinates = flip.axis.shiny(coordinates, x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), horizontal = T)
    coordinates$id = seurat$id

    img = seurat@images[[image]]@image
    img_grob = grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
    annotation = annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img))

    return(list(annotation, coordinates, img))
  }

  ai.filter = function(image = NULL, df = NULL, python_env = NULL, script = NULL, barcodeNum = NULL, thre = NULL, prefix = NULL) {
    writePNG(image, target = paste(getwd(), "STvis.AI.filter.genger.temp.png", sep = "/"))
    system(paste(python_env, script,
                 paste0("--img=", paste(getwd(), "STvis.AI.filter.genger.temp.png", sep = "/")),
                 paste0("--barcodeNum=", barcodeNum),
                 paste0("--thre=", thre), sep = " "))
    ai.filtered.pixels = read.csv(paste("~/STvis.AI.filter.genger.filterd.pixels.csv", sep = "/"), header = F)
    colnames(ai.filtered.pixels) = "barcode"
    ai.filtered.pixels$barcodeB = str_split(ai.filtered.pixels$barcode, "x", simplify = T)[ , 1]
    ai.filtered.pixels$barcodeA = str_split(ai.filtered.pixels$barcode, "x", simplify = T)[ , 2]
    ai.filtered.pixels$barcodeA = barcodeNum + 1 - as.numeric(ai.filtered.pixels$barcodeA)
    ai.filtered.pixels$barcode.use = paste0(ai.filtered.pixels$barcodeA, "x", ai.filtered.pixels$barcodeB)

    # add important feature!
    df[["ai.filter"]][df[["barcode"]] %in% paste0(prefix, ai.filtered.pixels$barcode.use)] = "filtered"

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

    if (length(Images(seurat)) > 1) {
      message(paste("Detect", length(Images(seurat)), "images!"))
      seurat.backup = seurat
      if (!is.null(image)) {
        seurat = subset(seurat.backup, cells = rownames(GetTissueCoordinates(seurat.backup, image = image)), slot = slot)
        seurat@meta.data = seurat@meta.data[ , colnames(seurat.backup@meta.data)]
        img = list(seurat@images[[image]]); names(img) = image
        seurat@images = img
      } else {
        warning("Please select your wanted image!")
        image = Images(seurat.backup)[1]
        seurat = subset(seurat.backup, cells = rownames(GetTissueCoordinates(seurat.backup, image = image)), slot = slot)
        seurat@meta.data = seurat@meta.data[ , colnames(seurat.backup@meta.data)]
        img = list(seurat@images[[image]]); names(img) = image
        seurat@images = img
      }
    }

    if (sum(str_detect(colnames(seurat), "x")) < dim(seurat)[2]) {
      stop("Your cell id must be formatted like 1x1, 1x2 or sample_1sx1, sample_1x2 ...")
    }

    sampleChoice = as.character(unique(seurat$orig.ident))
    shapeChoice = c(22, 21)
    featureChoice = colnames(seurat@meta.data)
    # for (i in featureChoice) {
    #   seurat@meta.data[is.na(seurat[[i]]), i] = "STvis"
    #   seurat@meta.data[is.null(seurat[[i]]), i] = "STvis"
    # }

    clean.barcodes = str_match(colnames(seurat), pattern = "\\d+x\\d+")[ , 1]
    prefix = gsub(clean.barcodes[1], "", colnames(seurat)[1])
    # add important feature!
    seurat$barcodeB = str_split(clean.barcodes, "x", simplify = T)[ , 1]
    seurat$barcodeA = str_split(clean.barcodes, "x", simplify = T)[ , 2]
    seurat$id = 1:dim(seurat)[2]

    updateSelectInput(session, inputId = "sampleInput", label = "Select sample", choices = sampleChoice, selected = sampleChoice[1])
    updateSelectInput(session, inputId = "shapeInput", label = "Select shape", choices = shapeChoice, selected = 22)
    updateSelectInput(session, inputId = "featureInput", label = "Select feature", choices = featureChoice, selected = "orig.ident")
    updateSelectInput(session, inputId = "subsetFeature", label = "Select feature to subset", choices = featureChoice, selected = "orig.ident")

    df = reactiveValues()
    df.backup = reactiveValues()
    for (i in c(featureChoice, "id")) {
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
      rv$ann = add_image(seurat)
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
      rv$ann[[2]] = flip.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), horizontal = T)
    })

    observeEvent(input$horizontal, {
      rv$ann[[2]] = flip.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), vertical = T)
    })

    observeEvent(input$angle.ok, {
      rv$ann[[2]] = rotate.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), angle = input$angle)
    })

    observeEvent(input$xMove.ok, {
      rv$ann[[2]] = move.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), x.num = input$xMove)
    })

    observeEvent(input$yMove.ok, {
      rv$ann[[2]] = move.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), y.num = input$yMove)
    })

    observeEvent(input$xShrink.ok, {
      rv$ann[[2]] = shrink.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), x.scale.factor = input$xShrink)
    })

    observeEvent(input$yShrink.ok, {
      rv$ann[[2]] = shrink.axis.shiny(rv$ann[[2]], x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), y.scale.factor = input$yShrink)
    })

    output$Plot1 = renderGirafe({
      withProgress(message = "Updating plot", value = 0,
                   {
                     df.tmp = as.data.frame(reactiveValuesToList(df))
                     x = girafe(ggobj = make.feature.plot.shiny(ann = rv$ann, anno.df = df.tmp, alpha = input$alphaValue, pt.size = input$spotSize*10, mode = input$filter.mode,
                                                                shape = as.integer(input$shapeInput), show.feature = ifelse(is.null(input$featureInput), "orig.ident", input$featureInput)),
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
      df[["ai.filter"]] = rep("exist", length(df[["id"]]))
      df.backup[["ai.filter"]] = rep("exist", length(df.backup[["id"]]))
      featureChoice = unique(c(featureChoice, "ai.filter"))
      updateSelectInput(session, inputId = "featureInput", label = "Select feature", choices = featureChoice, selected = "ai.filter")
      updateSelectInput(session, inputId = "subsetFeature", label = "Select feature to subset", choices = featureChoice, selected = "ai.filter")

      df[["ai.filter"]] = ai.filter(rv$ann[[3]], df, python_env = python_env, script = script, barcodeNum = ifelse(max(seurat$barcodeB) > 50, 96, 50), thre = input$thre, prefix = prefix)
      df.backup[["ai.filter"]] = ai.filter(rv$ann[[3]], df.backup, python_env = python_env, script = script, barcodeNum = ifelse(max(seurat$barcodeB) > 50, 96, 50), thre = input$thre, prefix = prefix)

      session$sendCustomMessage(type = "Plot1_set", message = character(0))
    })

    observeEvent(input$confirm, {
      shinyalert("Remember to click [quit] to save your confirmed labels!")
      ids.selected = as.numeric(input$Plot1_selected)
      if (is.factor(df[[input$featureInput]])) {
        levs = levels(df[[input$featureInput]])
        df[[input$featureInput]] = as.character(df[[input$featureInput]])
        df[[input$featureInput]][which(df$id %in% ids.selected)] = input$labelInput
        df[[input$featureInput]] = factor(df[[input$featureInput]], levels = unique(c(levs, input$labelInput)))

        levs = levels(df.backup[[input$featureInput]])
        df.backup[[input$featureInput]] = as.character(df.backup[[input$featureInput]])
        df.backup[[input$featureInput]][which(df.backup$id %in% ids.selected)] = input$labelInput
        df.backup[[input$featureInput]] = factor(df.backup[[input$featureInput]], levels = unique(c(levs, input$labelInput)))
      } else {
        df[[input$featureInput]][which(df$id %in% ids.selected)] = input$labelInput

        df.backup[[input$featureInput]][which(df.backup$id %in% ids.selected)] = input$labelInput
      }
      session$sendCustomMessage(type = "Plot1_set", message = character(0))
    })

    observe({

      if (input$stopApp > 0) {
        spots.df = df.backup$barcode
        spots.seurat = rownames(seurat@meta.data)

        featureChoice = featureChoice[!featureChoice %in% c("barcodeB", "barcodeA", "id")]

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
          seurat.backup@images[[image]]@coordinates[ , c("imagerow", "imagecol")] = rv$ann[[2]][ , c("x", "y")] %>%
            rotate.axis.shiny(x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), angle = 90) %>%
            flip.axis.shiny(x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), horizontal = T)
          for (i in colnames(seurat.backup@meta.data)) {
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

          stopApp(returnValue = seurat.backup)
        } else {
          seurat@images[[image]]@coordinates[ , c("imagerow", "imagecol")] = rv$ann[[2]][ , c("x", "y")] %>%
            rotate.axis.shiny(x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), angle = 90) %>%
            flip.axis.shiny(x = "x", y = "y", numBarcode = ifelse(max(seurat$barcodeB) > 50, 96, 50), horizontal = T)
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

