
##' Handle external browser options for the shiny.
##'
##'
##' @import shiny
##' @export
##' @param openinBrowser a logical indicating whether the shiny app should be opened in the browser or not
browserManager <- function(openinBrowser){
  if (openinBrowser){
    ##Use get to circumvent the global binding errors.
    options(shiny.launch.browser = get(".rs.invokeShinyWindowExternal", envir = as.environment("tools:rstudio")))
  }
  else{
    options(shiny.launch.browser = get(".rs.invokeShinyWindowViewer", envir = as.environment("tools:rstudio")))
  }
}

loaders <- c("Serum Metabolomics",
             "Gut Microbiome",
             "Iglu CGM",
             "Liver Ultrasound",
             "ABI",
             "Sleep",
             "Hormonal Status",
             "DEXA")

get_basepath <- function(loader, use_clumped){

  if (loader == "Serum Metabolomics"){
    if (use_clumped == "Yes"){
      basepath <- "/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/metab/gwas_results_clumped"
    }
    else{
      basepath <- "/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/metab/gwas_results_metab"
    }
  }
  else if (loader == "Gut Microbiome"){
    if (use_clumped == "Yes"){
      basepath <- "/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/microbiome/gwas_results_clumped"
    }
    else{
      basepath <- "/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/microbiome/gwas_results_mb"
    }
  }
  else{
    if (use_clumped== "Yes"){
      basepath <- "/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results_clumped"
    }
    else{
      basepath <- "/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results"
    }
  }
  basepath
}

get_available_gwases <- function(loader, use_clumped){
  basepath <- get_basepath(loader, use_clumped)
  thefiles<- list.files(basepath)
  thefiles <- thefiles[thefiles != "batch0.log"]
  thefiles
}

##from https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
ggplot_manhattan <- function(gwas_data, theTitle = "Manhattan Plot", ymin = 4){
  gwas_data <- gwas_data[, colnames(gwas_data) %in% c("BP", "P", "CHR", "SNP")] ##fix empty column name errors
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

  sig <- 5e-8
  manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P),
                                    color = as.factor(CHR), size = -log10(P))) +
    geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
    geom_point(alpha = 0.75) +
    scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
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

##' Run the McMasterPandemic Shiny
##'
##' run_shiny() is an example of a single-file Shiny app in that it defines a UI object and a server method to handle that object. A benefit is that Roxygen import tags only need to be called once to become available to the entire shiny, and it's easy to run the shiny as well. After both the ui object and server function are defined, browserManager is called, which sets the viewing environment that the Shiny runs in. Anywhere with a tag$html call is either a CSS or HTML tag to change certain visual aspects of the shiny.  In addition, many ui elements are rendered in the server function and passed to ui with renderUI. RenderHTML is another wrapper for this as well. This lets us make UI elements which depend on input from the UI itself.
##'
##' @param useBrowser A logical indicating whether to run the Shiny in the default browser or in a seperate window generated by RStudio.
##' @param usingOnline A logical indicating whether the app is being run on shinyapps.io or not.
##' @import shiny
##' @return NULL
##' @export
run_shiny <- function(useBrowser = TRUE, usingOnline = FALSE) {
  ui <- fluidPage(theme = shinythemes::shinytheme("flatly"),
                  ##Set the title panel to be Heritage Maroon.
                  h1(id = "heading", "Eran Segal 10K Project Interactive GWAS Results Interface"),
                  ##Colour the top and bottom of the page appropriately.
                  tags$style(HTML("#heading {background-color: #0078a4; color: white !important;}")),
                  tags$style(HTML("#sourcelink {background-color: #0078a4; color: white !important;}")),
                  #Change the colours of text on the tab selectors to be blue.
                  tags$style(HTML("
                                  .tabbable > .nav > li > a[data-value='summarystatspanel'] {color: black}
                                  .tabbable > .nav > li > a[data-value='gencorr'] {color: black}
                                  ")),
                  shinyWidgets::setBackgroundColor(color = "#e6ebed"),
                  sidebarLayout(
                    sidebarPanel(id = "sidebar", width = 4,
                                 selectInput("do_clumped",
                                             label = "Use clumped results for LD",
                                             choices = c("Yes", "No")),
                                 ## Only show the selector to input parameters if that's selected.
                                   tabsetPanel(
                                     id = "tabs",
                                     selected = "summarystatspanel",
                                     tabPanel(
                                       title = "Download ldscore report",
                                       value = "gencorr",
                                       selectInput("domain1",
                                                   label = "Phenotype 1 domain:",
                                                   choices = loaders),
                                       uiOutput("pheno1"),
                                       selectInput("domain2",
                                                   label = "Phenotype 2 domain:",
                                                   choices = loaders),
                                       uiOutput("pheno2"),
                                       uiOutput("ldDownload")
                                     ),
                                     tabPanel(
                                       title = "Download GWAS Summary Statistics",
                                       value = "summarystatspanel",
                                       selectInput("domain",
                                                   label = "Data domain:",
                                                   choices = loaders),
                                       uiOutput("pheno"),
                                       uiOutput("gwasDownload")
                                       ))),
                    mainPanel(
                      fluidRow(
                        uiOutput("plotColumn"),
                        br())
                    )
                  ),
                  uiOutput("sourcelink")
                  )

  #Everything else.
  server <- function(input, output, session){
    output$plot <- renderPlot({
      qqman::manhattan(read_clumped(paste0(get_basepath(input$domain, input$do_clumped), "/", input$pheno)))
      })
    output$plotColumn <- renderUI({
      plotOutput({"plot"})
    })
    output$pheno <- renderUI({
      selectInput("pheno",
                  label = "Phenotype:",
                  choices = get_available_gwases(input$domain, input$do_clumped))
    })
    output$gwasDownload <- renderUI({
      downloadButton("downloadData", "Download GWAS Summary Statistics", class = "dbutton")
    })
    output$pheno1 <- renderUI({
      selectInput("pheno1",
                  label = "Phenotype:",
                  choices = get_available_gwases(input$domain1, input$do_clumped))
    })
    output$pheno2 <- renderUI({
      selectInput("pheno2",
                  label = "Phenotype:",
                  choices = get_available_gwases(input$domain2, input$do_clumped))
    })
    output$ldDownload <- renderUI({
      downloadButton("downloadData", "Download ldscore report", class = "dbutton")
    })


    ##Handle downloads for the sample template csv file.
    ##The file is an empty version of ICU1.csv so it scales as more parameters are added.
    output$downloadData <- downloadHandler(
      filename = function(){return("gwas_results")},
      content = function(file) {
        basicTemplate <- read_in_csv()
        ##Default the values column.
        basicTemplate$Value <- read_params("PHAC.csv")
        write.csv(basicTemplate, file, row.names = FALSE)
      }
    )
      output$sourcelink <- renderUI({tagList(
      "Source code available at ",
      a("https://github.com/repo-link", href = "https://github.com/repo-link"),
      "                                GWAS Repository"
    )})
  }
  if (usingOnline == FALSE){
    ##Set the viewing options first.
    ##Run the shiny app. the default value of launch.browser looks for the option set by browserManager.
    shiny::runApp(appDir = shinyApp("ui" = ui, "server" = server), launch.browser = getOption("shiny.launch.browser", interactive()))
    browserManager(useBrowser)
  }
  else{
    shinyApp(ui, server)
    }
}
