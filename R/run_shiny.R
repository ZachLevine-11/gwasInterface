

##' @export
##' @import fastman
shiny_manhattan <- function(x, phenoName){
  ##ylims: -log10(5*10**-n) ~ n-1, i.e for 4 put 3
  fastman::fastman(x, chr = "X.CHROM", bp ="POS", cex.axis = 1, annotatePval =  5*10**(-8)/727, snp = "ID", cex.text = 1.3, genomewideline = -log10(5*10**(-8)/727), suggestiveline = -log10(5e-8), ylim = c(min(2.5, min(-log10(x$P))), max(-log10(x$P))), main = phenoName)

}

options_mb <-   read.csv(system.file("lists",
                                     "mb_options.csv",
                                     package = "gwasInterface"))
list_kingdom <- function(){
  c(unique(options_mb$kingdom))

}
list_phylum <- function(kingdom){
  print(kingdom)
  unique(options_mb[options_mb$kingdom == kingdom, "phylum"])
}
list_class <- function(kingdom, phylum){
  temp <- options_mb[options_mb$kingdom == kingdom,]
  unique(temp[temp$phylum == phylum, "class"])

}
list_order <- function(kingdom, phylum, class){
  temp <- options_mb[options_mb$kingdom == kingdom,]
  temp2 <- temp[temp$phylum == phylum,]
  unique(temp2[temp2$class == "class","order"])
}
list_family <- functiony <- function(kingdom,phylum, class, order){
  temp <- options_mb[options_mb$kingdom == kingdom,]
  temp2 <- temp[temp$phylum == phylum,]
  temp3 <- temp2[temp2$class == class,]
  unique(temp[temp$order == order, "family"])
}
list_genus <- function(kingdom, phylum, class, order, family){
  temp <- options_mb[options_mb$kingdom == kingdom,]
  temp2 <- temp[temp$phylum == phylum,]
  temp3 <- temp2[temp2$class == class,]
  temp4 <- temp3[temp$order == order,]
  unique(temp4[temp4$family == family, "genus"])

}
list_species <- function(kingdom, phylum, class,order, family, genus){
  temp <- options_mb[options_mb$kingdom == kingdom,]
  temp2 <- temp[temp$phylum == phylum,]
  temp3 <- temp2[temp2$class == class,]
  temp4 <- temp3[temp$order == order,]
  temp5 <- temp4[temp4$family == family,]
  unique(temp5[temp5$genus == genus, "species"])
}

find_fname <- function(kingdom, phylum, class,order, family, genus, species){
  temp <- options_mb[options_mb$kingdom == kingdom,]
  temp2 <- temp[temp$phylum == phylum,]
  temp3 <- temp2[temp2$class == class,]
  temp4 <- temp3[temp$order == order,]
  temp5 <- temp4[temp4$family == family,]
  temp6<- (temp5[temp5$genus == genus,])
  temp6[[temp5$species == species,"X"]]
  }


list_features <- function(domain){
  read.csv(system.file("lists",
                       paste0(domain ,"Loader", ".csv"),
                       package = "gwasInterface"))[,2]
}

make_paper_figures <- function(){
  qu <- read.csv("/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results/batch0.1st_Qu..glm.linear", sep = "\t")
 # dx <- read.csv("/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results/batch0.body_head_bmd.glm.linear", sep = "\t")
  #jpeg("~/Desktop/dx.jpeg", width = 960)
  #fastman(dx, chr = "X.CHROM", bp ="POS", ylim = c(0, 11), cex.axis = 1, annotatePval =  5*10**(-8)/727, snp = "ID", cex = 1, genomewideline = -log10(5*10**(-8)/727), suggestiveline = -log10(5e-8))
  #dev.off()
  jpeg("~/Desktop/qu.jpeg", width = 960)
  fastman(qu, chr = "X.CHROM", bp ="POS", cex.axis = 1, annotatePval =  5*10**(-8)/727, snp = "ID", cex.text = 1.3, genomewideline = -log10(5*10**(-8)/727), suggestiveline = -log10(5e-8))
  dev.off()
}

make_cnv_manplot <- function(){
  library(qqman)
  res <- read.csv("/net/mraid08/export/jasmine/zach/cnv/ea1c_test_formatted.csv")
  res$CHR <- as.integer(res$CHR)
  res_add_to_each_chr_start <- sapply(unique(res$CHR), function(chr){res[res$CHR == toString(as.numeric(chr) - 1), "END"][length(res[res$CHR == toString(as.numeric(chr) - 1),])]})
  res_add_to_each_chr_start[is.na(res_add_to_each_chr_start)] <- 0
  res_add_to_each_chr_start <- cumsum(res_add_to_each_chr_start)
  for (chr in names(res_add_to_each_chr_start)){
    res[res$CHR == toString(chr),"START"] <- res[res$CHR == toString(chr),"START"] + res_add_to_each_chr_start[chr]
  }
  res$BP <- (res$START + res$END)/2
  res[res$CHR == "X", "CHR"] <- 23
  res[res$CHR == "Y", "CHR"] <- 24
  res[res$CHR == "MT", "CHR"] <- 25
  res$CHR <- as.numeric(res$CHR)
  jpeg("~/Desktop/ea1c_cnv.jpeg", width = 960)
  manhattan(res[!is.na(res$P),])
  dev.off()
  }

ea1c_plot <- function(){
  library(fastman)
  ea1c <- read.csv("/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results/batch0.ea1c.glm.linear", sep = "\t")
  jpeg("~/Desktop/ea1c_snps.jpeg", width = 960)
  fastman(ea1c, chr = "X.CHROM", bp ="POS", cex.axis = 1, annotatePval =  5*10**(-8)/727, snp = "ID", cex.text = 1.3, genomewideline = -log10(5*10**(-8)/727), suggestiveline = -log10(5e-8))
  dev.off()
  }

##from https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
ggplot_manhattan <- function(gwas_data, theTitle = "Manhattan Plot", ymin = 4){
  colnames(gwas_data) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1", "AX", "TEST", "N", "BETA", "SE", "T_STAT", "P")
  data_cum <- gwas_data %>%
    group_by(CHR) %>%
    summarise(max_bp = max(BP)) %>%
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
    select(CHR, bp_add)
  gwas_data <- gwas_data %>%
    inner_join(data_cum, by = "CHR") %>%
    mutate(bp_cum = BP + bp_add)
  axis_set <- gwas_data %>%
    group_by(CHR) %>%
    summarize(center = mean(bp_cum))

  ylim <- gwas_data %>%
    filter(P == min(P)) %>%
    mutate(ylim = abs(floor(log10(P))) + 2) %>%
    pull(ylim)

  manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P),
                                    color = as.factor(CHR), size = -log10(P))) +
    geom_point(alpha = 0.75) +
    scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center, guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(expand = c(0,0), limits = c(ymin, ylim)) +
    scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$CHR)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL,
         #https://sites.google.com/view/stuck-in-the-shallow-end/home/generate-manhattan-plots-with-ggplot2-and-ggrastr
         title = theTitle) +
    ylab(expression(paste(-log[10],"(", italic(P), ")"))) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    )
  manhplot <- manhplot  + theme(title = element_text(size = 0), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 17), axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 15)) + labs(x = "CHR")
  manhplot
}

read_clumped <- function(fname){
  read_in <- do.call(rbind, lapply(strsplit(readLines(fname), "\\s+|\\t+|\\s+\\t+|\\t+\\s+"), function(x){as.data.frame(t(x))}))
  read_in <- read_in[, 1:6]
  colnames(read_in) <- read_in[1,]
  read_in <- read_in[-1,] ##Drop the column names as the first entry
  read_in$CHR <- as.numeric(read_in$CHR)
  read_in$"F" <- as.numeric(read_in$"F")
  read_in$BP <- as.numeric(read_in$BP)
  read_in$P <- as.numeric(read_in$P)
  read_in
}

format_gwas <- function(thephenoname){
  thegwas <- readRDS(system.file(paste0("full_results_rds/", thephenoname, ".Rds"),
                                 package = "gwasInterface"))
  colnames(thegwas) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1", "AX", "TEST", "N", "BETA", "SE", "T_STAT", "P")
  thegwas <- thegwas[order(thegwas$P),]
  data.table(thegwas)
}

domains_list <- c("CGM", "Ultrasound", "HormonalStatus", "DEXA", "RetinaScan", "Metabolomics", "Microbiome")

##' Run the GWASInterface Shiny
##'
##' run_shiny() is an example of a single-file Shiny app in that it defines a UI object and a server method to handle that object. A benefit is that Roxygen import tags only need to be called once to become available to the entire shiny, and it's easy to run the shiny as well. After both the ui object and server function are defined, browserManager is called, which sets the viewing environment that the Shiny runs in. Anywhere with a tag$html call is either a CSS or HTML tag to change certain visual aspects of the shiny.  In addition, many ui elements are rendered in the server function and passed to ui with renderUI. RenderHTML is another wrapper for this as well. This lets us make UI elements which depend on input from the UI itself.
##'
##' @import shiny
##' @import ggplot2
##' @import shinythemes
##' @import dplyr
##' @import tibble
##' @import data.table
##' @return NULL
##' @export
run_shiny <- function(useBrowser = TRUE, usingOnline = FALSE) {
  ui <- fluidPage(theme = shinythemes::shinytheme("flatly"),
                  h1(id = "heading", "Eran Segal Human Phenotype Project Interactive GWAS Results Interface"),
                  ##Colour the top and bottom of the page appropriately.
                  tags$style(HTML("#heading {background-color: #0078a4; color: white !important;}")),
                  tags$style(HTML("#sourcelink {background-color: #0078a4; color: white !important;}")),
                  #Change the colours of text on the tab selectors to be blue.
                  tags$style(HTML("
                                  .tabbable > .nav > li > a[data-value='summarystatspanel'] {color: black}
                                  .tabbable > .nav > li > a[data-value='gencorr'] {color: black}
                                  ")),
#                  shinyWidgets::setBackgroundColor(color = "#e6ebed"),
                  sidebarLayout(
                    sidebarPanel(id = "sidebar", width = 4,
                                 ## Only show the selector to input parameters if that's selected.
                                   tabsetPanel(
                                     id = "tabs",
                                     selected = "summarystatspanel",
                                     tabPanel(
                                       uiOutput("pheno1"),
                                     ),
                                     tabPanel(
                                       title = "Feature Selection",
                                       value = "summarystatspanel",
                                       selectInput("domain",
                                                   label= "Feature Domain",
                                                   choices = domains_list),
                                       h4(id = "t0", "PRS Associations Heatmap (Corrected P)"),
                                       plotOutput("heatmap"),
                                       uiOutput("featureSelect"),
                                       uiOutput("aux_featureSelect"),
                                       uiOutput("gwasDownload")
                                       ))),
                    mainPanel(
                        h3(id = "t1", "Manhattan Plot of GWAS Results"),
                        plotOutput("plot", width  = "100%"),
                        ##Colour the top and bottom of the page appropriately.
                        h3(id = "t2", "Top Multi-SNP Significant PRS Associations"),
                        dataTableOutput("prs_assoc"),
                        h3(id = "t3", "Interactive Table of GWAS Results"),
                        dataTableOutput("table")
                    )
                  ),
                  uiOutput("sourcelink")
                  )

  #Everything else.
  server <- function(input, output, session){
    prs_table_loaded <- readRDS(system.file("prs_assoc_corrected_clinical.Rds",
                                             package = "gwasInterface"))
    prs_name_dict_loaded <- readRDS(system.file("prs_name_dict.Rds",
                                            package = "gwasInterface"))
    observeEvent(input[["pheno"]], {
      output$plot <- renderPlot({
        ##Force reactive loading based on these values
        if (input$domain == "Microbiome"){
          actual_pheno_fname = find_fname <-find_fname(input$kingdom, input$phylum, input$class,input$order, input$family, input$genus,input$species)
          shiny_manhattan(readRDS(system.file(paste0("full_results_rds/",actual_pheno_fname, ".Rds"),
                                              package = "gwasInterface")), input$pheno)
        }
        else{
          shiny_manhattan(readRDS(system.file(paste0("full_results_rds/", input$pheno, ".Rds"),
                                              package = "gwasInterface")), input$pheno)
        }
        })

      output$table <- renderDataTable({
        thegwas <- format_gwas(input$pheno)
        })
      output$prs_assoc <- renderDataTable({
        justpheno <- as.numeric(prs_table_loaded[prs_table_loaded$Phenotype == input$pheno,]) ##first two entries are NA corresponding to the index of Loader, Phenotype
        sig_ids <- justpheno < 0.05
        sig_ids[c(1,2)] <- c(FALSE, FALSE) ##set FALSE (instead of NA) to the first two columns because they are just the index
        justpheno_sig_P <- justpheno[sig_ids]
        justpheno_sig_prs_names <- colnames(prs_table_loaded)[sig_ids]
        ##use a tibble so we can have spaces in the column names
        res_table <- tibble("Phenotype" = rep(input$pheno, length(justpheno_sig_P)), "Corrected P Value" = justpheno_sig_P, "Polygenic Risk Score (PRS)" = justpheno_sig_prs_names)
        res_table["Corresponding UK Biobank Phenotype Code for PRS"] <- sapply(res_table[["Polygenic Risk Score (PRS)"]], function(thename){prs_name_dict_loaded[prs_name_dict_loaded$name == thename,]$code})
        res_table <- data.table(res_table[order(res_table["Corrected P Value"]), c("Phenotype", "Polygenic Risk Score (PRS)", "Corresponding UK Biobank Phenotype Code for PRS", "Corrected P Value")])
      })
    })

    observeEvent(input[["domain"]], {
        output$heatmap <- renderImage({
          list(src=system.file(paste0("prs_heatmap_imgs/", input$domain, "Loader", ".png"),
                               package = "gwasInterface"),
               width = "80%",
               height = "100%")
        }, deleteFile = FALSE)
    })

    output$gwasDownload <- renderUI({
      downloadButton("downloadData", "Download GWAS Summary Statistics", class = "dbutton")
    })
    output$featureSelect <- renderUI({
      if (input$domain != "Microbiome"){

      }

    })
    output$aux_featureSelect <- renderUI({
        if (input$domain == "Microbiome"){
        selectInput("kingdom",
                    label = "Kingdom:",
                    choices = list_kingdom(),
                    selected = "Bacteria")
          selectInput("phylum",
                      label = "Phylum:",
                      choices = list_phylum(input$kingdom))


        selectInput("class",
                    label = "Class:",
                    choices = list_class(input$kingdom, input$phylum))
        selectInput("order",
                    label = "Order:",
                    choices = list_order(input$kingdom, input$phylum, input$class))
        selectInput("family",
                    label = "Family:",
                    choices = list_family(input$kingdom, input$phylum, input$class, input$order))
        selectInput("genus",
                    label = "Genus:",
                    choices = list_genus(input$kingdom, input$phylum, input$class, input$order, input$family))
        selectInput("species",
                    label = "Species:",
                    choices = list_species(input$kingdom, input$phylum, input$class, input$order, input$family, input$genus))
        }
      })

    ##Handle downloads for the sample template csv file.
    ##The file is an empty version of ICU1.csv so it scales as more parameters are added.
    output$downloadData <- downloadHandler(
      filename = function(){return(paste0(input$pheno, ".csv"))},
      content = function(file) {
        write.csv(format_gwas(input$pheno), file, row.names = FALSE)
      }
    )
      output$sourcelink <- renderUI({tagList(
      "Source code available at ",
      a("https://github.com/ZachLevine-11/gwasInterface", href = "https://github.com/ZachLevine-11/gwasInterface"),
      "                                Based on the paper: Paper title"
    )})
  }
    shinyApp(ui, server)
}
