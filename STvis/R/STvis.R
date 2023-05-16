library(magick)
library(shiny)
library(shinydashboard)
library(png)
library(tools)
library(DT)

library(Seurat)
library(dplyr)
library(tibble)
library(data.table)
library(stringr)

library(ggplot2)
library(viridis)
library(grid)

shiny_st = function(seurat, assay = "SCT", slot = "data") {
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

  make.feature.plot.shiny = function(ann = NULL, anno.df = NULL, alpha = 0.8, pt.size = 0.1, shape = 22, show.feature = NULL) {
    annotation = ann[[1]]
    coordinates = ann[[2]]
    img = ann[[3]]

    coordinates = coordinates[paste0(anno.df$barcodeB, "x", anno.df$barcodeA), ]
    coordinates$feature = anno.df[[show.feature]]

    if (is.numeric(coordinates$feature)) {
      cols = viridis(100, option = "D")
      p = ggplot(coordinates, aes(x = x, y = y, data_id = id)) +
        annotation +
        ggiraph::geom_point_interactive(aes(fill = feature, alpha = feature), size = pt.size, shape = shape, stroke = 0) +
        scale_fill_gradientn(colors = cols) +
        scale_alpha(range = c(alpha, 1)) +
        ylim(nrow(img), 0) + xlim(0, ncol(img)) +
        theme_void() + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = T, clip = "on") +
        theme(legend.position = "top") +
        guides(alpha = "none") +
        theme(aspect.ratio = 1) +
        labs(fill = show.feature)
    } else {
      cols = c("#5A5156FF", "#E4E1E3FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", "#2ED9FFFF", "#DEA0FDFF", "#AA0DFEFF",
               "#F8A19FFF", "#325A9BFF", "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", "#782AB6FF",
               "#AAF400FF", "#BDCDFFFF", "#822E1CFF", "#B5EFB5FF", "#7ED7D1FF", "#1C7F93FF", "#D85FF7FF", "#683B79FF", "#66B0FFFF", "#3B00FBFF")
      p = ggplot(coordinates, aes(x = x, y = y, data_id = id)) +
        annotation +
        ggiraph::geom_point_interactive(aes(fill = feature), size = pt.size, shape = shape, stroke = 0, alpha = alpha) +
        scale_fill_manual(values = cols) +
        ylim(nrow(img), 0) + xlim(0, ncol(img)) +
        theme_void() + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = T, clip = "on") +
        theme(legend.position = "top") +
        guides(alpha = "none") +
        theme(aspect.ratio = 1) +
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

  ui = dashboardPage(dashboardHeader(title = tags$div(style = "white-space: pre-wrap; word-wrap: break-word; line-height: 1.2;",
                                                      "ST Visualization Tool (STvis v1.0.0)\n---please watch the [Instructions] first---"), titleWidth = 450),
                     dashboardSidebar(width = 600,
                                      actionButton(inputId = "info", label = "Instructions"), shiny::hr(),

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
                                          actionButton(inputId = "recover.ok", label = "Back to all celltypes!)", width = "50%"))), shiny::hr(),

                                      fluidRow(
                                        column(width = 6,
                                               sliderInput(inputId = "alphaValue", label = "Spot.alpha [0-1]", min = 0, max = 1, value = 0.8, step = 0.1)),
                                        column(width = 6,
                                               sliderInput(inputId = "spotSize", label = "Spot.size [0-1]", min = 0, max = 1, value = 0.1, step = 0.1))), shiny::hr(),

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
                     dashboardBody(ggiraph::girafeOutput("Plot1", width = "100%", height = paste0(1080, "px"))),
                     tags$head(
                       tags$style(HTML(".sidebar {width: 600px;}
                                 .my-row .form-group {margin-top: 0; margin-bottom: 0;}
                                 .content-wrapper {margin-left: 0;}
                                 .input-group {margin-right: 0;}")))
  )

  server = function(input, output, session) {
    options(shiny.maxRequestSize = 100*1024^2)

    if (sum(str_detect(colnames(seurat), "x")) < dim(seurat)[2]) {
      stop("Your cell name must be formatted like 1x1, 1x2 ...")
    }

    sampleChoice = unique(seurat$orig.ident)
    shapeChoice = c(22, 21)
    featureChoice = colnames(seurat@meta.data)
    # add important feature!
    seurat$barcodeB = str_split(colnames(seurat), "x", simplify = T)[ , 1]
    seurat$barcodeA = str_split(colnames(seurat), "x", simplify = T)[ , 2]
    seurat$id = 1:dim(seurat)[2]

    updateSelectInput(session, inputId = "sampleInput", label = "Select sample", choices = sampleChoice, selected = sampleChoice[1])
    updateSelectInput(session, inputId = "shapeInput", label = "Select shape", choices = shapeChoice, selected = 22)
    updateSelectInput(session, inputId = "featureInput", label = "Select feature", choices = featureChoice, selected = "orig.ident")
    updateSelectInput(session, inputId = "subsetFeature", label = "Select feature to subset", choices = featureChoice, selected = "orig.ident")

    df = reactiveValues()
    df.backup = reactiveValues()
    for (i in featureChoice) {
      df[[i]] = seurat@meta.data[[i]]
      df.backup[[i]] = seurat@meta.data[[i]]
    }

    rv = reactiveValues(sNr = "1", ann = NULL)
    observeEvent(input$sampleInput, {
      rv$sNr = sampleChoice[1]
      rv$ann = add_image(seurat)
    })

    observeEvent(input$gene.ok, {
      mat.df = FetchData(srat.merge, vars = input$geneInput, slot = slot) %>% rownames_to_column(var = "barcode")
      df.backup[[input$geneInput]] = mat.df[[input$geneInput]][mat.df$barcode %in% paste0(df.backup$barcodeB, "x", df.backup$barcodeA)]
      df[[input$geneInput]] = mat.df[[input$geneInput]][mat.df$barcode %in% paste0(df$barcodeB, "x", df$barcodeA)]
      updateSelectInput(session, inputId = "featureInput", label = "Select feature", choices = names(df), selected = input$geneInput)
      updateSelectInput(session, inputId = "subsetFeature", label = "Select feature to subset", choices = names(df), selected = "orig.ident")
    })

    observeEvent(input$subset.ok, {
      subset.features = strsplit(input$subsetString, split = ",")[[1]]

      if (all(subset.features %in% df[[input$subsetFeature]])) {
        for (i in names(df)[names(df) != input$subsetFeature]) {
          df[[i]] = df[[i]][df[[input$subsetFeature]] %in% subset.features]
        }
        df[[input$subsetFeature]] = df[[input$subsetFeature]][df[[input$subsetFeature]] %in% subset.features]
        if (is.factor(df[[input$subsetFeature]])) {
          df[[input$subsetFeature]] = droplevels(df[[input$subsetFeature]])
        }
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

    output$Plot1 = ggiraph::renderGirafe({
      withProgress(message = "Updating plot", value = 0,
                   {
                     df.tmp = as.data.frame(reactiveValuesToList(df))
                     x = ggiraph::girafe(ggobj = make.feature.plot.shiny(ann = rv$ann, anno.df = df.tmp, alpha = input$alphaValue, pt.size = input$spotSize*10,
                                                                         shape = as.integer(input$shapeInput), show.feature = ifelse(is.null(input$featureInput), "orig.ident", input$featureInput)),
                                         width_svg = 12, height_svg = 10)
                     x = ggiraph::girafe_options(x, ggiraph::opts_zoom(max = 6),
                                                 ggiraph::opts_selection(type = "multiple", css = "fill:transparent;stroke:transparent;opacity:0.7;"))
                     x
                   })
    })

    observeEvent(input$confirm, {
      ids.selected = as.numeric(input$Plot1_selected)
      if (is.factor(df[[input$featureInput]])) {
        levs = levels(df[[input$featureInput]])
        df[[input$featureInput]] = as.character(df[[input$featureInput]])
        df[[input$featureInput]][which(df$id %in% ids.selected)] = input$labelInput
        df[[input$featureInput]] = factor(df[[input$featureInput]], levels = c(levs, input$labelInput))

        levs = levels(df.backup[[input$featureInput]])
        df.backup[[input$featureInput]] = as.character(df.backup[[input$featureInput]])
        df.backup[[input$featureInput]][which(df.backup$id %in% ids.selected)] = input$labelInput
        df.backup[[input$featureInput]] = factor(df.backup[[input$featureInput]], levels = c(levs, input$labelInput))
      } else {
        df[[input$featureInput]][which(df$id %in% ids.selected)] = input$labelInput

        df.backup[[input$featureInput]][which(df.backup$id %in% ids.selected)] = input$labelInput
      }
      session$sendCustomMessage(type = "Plot1_set", message = character(0))
    })

    observe({
      if (input$stopApp > 0) {
        print("Stopped")
        spots.df = paste0(df$barcodeB, "x", df$barcodeA)
        spots.seurat = paste0(seurat$barcodeB, "x", seurat$barcodeA)
        for (i in featureChoice) {
          if (is.factor(seurat@meta.data[[i]])) {
            levs = levels(seurat@meta.data[[i]])
            seurat@meta.data[[i]] = as.character(seurat@meta.data[[i]])
            seurat@meta.data[[i]][spots.seurat %in% spots.df] = as.character(df[[i]])
            seurat@meta.data[[i]] = factor(seurat@meta.data[[i]], levels = c(levs, levels(df[[i]])[!levels(df[[i]]) %in% levs]))
          } else {
            seurat@meta.data[[i]][spots.seurat %in% spots.df] = df[[i]]
          }
        }
        stopApp(returnValue = seurat)
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

  runApp(list(ui = ui, server = server), launch.browser = getOption("shiny.launch.browser", interactive()))
}

