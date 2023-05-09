# Import libraries
library(rio)
library(shiny)
library(DT)
library(colourpicker)
library(gplots)
library(tidyverse)
library(ggbeeswarm)


# Define UI for app
ui <- fluidPage(
  
  # Application title
  titlePanel("mRNA Expression Analysis on Huntington's Disease vs Normal"),
  tags$p(" -- A Rshiny App by Yuxiang Zhou for Final Project for BF591 Spring 2023"),
  
  # Main tab panel for shiny app
  tabsetPanel(
    # Samples main-tab
    tabPanel("Samples",
             # a title with short description of this page's function
             tags$h3("Explore the metadata and the samples"),
             tags$div(
               tags$h5("Explore, summarize, search, and plot sample metadata")
             ),
             sidebarLayout( # Sidebar layout for samples tab
               sidebarPanel( # Sidebar panel for samples tab
                 fileInput(inputId = "samples", # Input for uploading samples file
                           label = "Please Upload Samples File", 
                           multiple = FALSE, 
                           accept = ".csv",
                           buttonLabel = "Browse Sample File",
                           placeholder = "data/metadata.csv"
                 )
               ),
               mainPanel( # Main panel for samples tab
                 tabsetPanel( # Sub-tabs panel for samples main-tab
                   tabPanel("Summary", # Summary sub-tab
                            tableOutput("sample_summary"),
                            htmlOutput("sample_info")
                   ),
                   tabPanel("Tables", # Tables sub-tab for samples main-tab
                            sidebarLayout(
                              sidebarPanel(
                                checkboxGroupInput(inputId = "sample_display_col", # Checkbox selection of columns to display
                                                   label = "Select Columns to Display:",
                                                   choices = list("Diseases Diagnosis" = "diagnosis:ch1", 
                                                                  "RIN Score" = "rin:ch1", 
                                                                  "PMI" = "pmi:ch1", 
                                                                  "mRNA Reads" = "mrna-seq reads:ch1", 
                                                                  "Age of Death" = "age of death:ch1", 
                                                                  "Age of Onset" = "age of onset:ch1", 
                                                                  "Duration of Disease" = "duration:ch1", 
                                                                  "CAG Repeat Length" = "cag:ch1", 
                                                                  "Vonsattel Grade" = "vonsattel grade:ch1", 
                                                                  "H-V Striatal Score" = "h-v striatal score:ch1", 
                                                                  "H-V Cortical Score" = "h-v cortical score:ch1"),
                                                   selected = c("diagnosis:ch1", "pmi:ch1", "age of death:ch1", "rin:ch1", "mrna-seq reads:ch1", "age of onset:ch1", "duration:ch1", "cag:ch1", "vonsattel grade:ch1", "h-v striatal score:ch1", "h-v cortical score:ch1"))
                                ,
                                actionButton(
                                  inputId = "sample_simple_table_display_button",
                                  label = "Display Selected Table"
                                ),
                                tags$hr(),
                                actionButton(
                                  inputId = "sample_table_display_button",
                                  label = "Display Full Table"
                                )
                              ),
                              mainPanel(
                                div(
                                  style = "font-family: Arial, sans-serif; font-size: 14px; line-height: 1.5;",
                                  tags$h5("A shortened table of the metadata is displayed below for simplicity."), 
                                  tags$h5("You can select the columns you would like to display."), 
                                  tags$h5("To display the full table, click the ", 
                                          tags$span("Display Full Table", style = "font-weight: bold;"),
                                          " button at left.")
                                ),
                                uiOutput("sample_table_display")
                              ),
                            )
                   ),
                   tabPanel("Plots", # Plots sub-tab for samples main-tab
                            sidebarLayout( # Sidebar layout for plots sub-tab
                              sidebarPanel( # Sidebar panel for plots sub-tab
                                fluidRow(
                                  column(
                                    width = 12,
                                    selectInput(
                                      inputId = "sample_plot_dropdown",
                                      label = "Select What to Draw:",
                                      choices = c("PMI", "RIN Score", "mRNA Reads", "Age of Death")
                                    ),
                                    tags$div(
                                      style = "margin-bottom: 5px;",
                                      actionButton(
                                        inputId = "plot_histogram",
                                        label = "Histogram"
                                      )
                                    ),
                                    tags$div(
                                      style = "margin-bottom: 5px;",
                                      actionButton(
                                        inputId = "plot_density",
                                        label = "Density Plot"
                                      )
                                    )
                                  )
                                )
                              ),
                              mainPanel( # Main panel for plots sub-tab
                                uiOutput("sample_plot_display")
                              )
                            )
                   )
                 )
               )
             )
    ),
    
    
    # Counts tab
    tabPanel("Counts",
            tags$h3("Explore the counts"),
            tags$div(
              tags$h5("Explore, summarize, search, and plot normalized counts")
            ),
             sidebarLayout( # Sidebar layout for counts main-tab
               sidebarPanel( # Sidebar panel for counts main-tab
                 fileInput(inputId = "counts", # Input for uploading counts file
                           label = "Please Upload Counts File", 
                           multiple = FALSE, 
                           accept = ".csv",
                           buttonLabel = "Browse Normalized Counts File",
                           placeholder = "data/norm_counts.csv"
                 ),
                 sliderInput(inputId = "count_percentile_slider", 
                             label = "Genes With at Least `X` Percentile of Variance:", 
                             min = 1, 
                             max = 100, 
                             value = 80
                 ),
                 sliderInput(inputId = "count_nonzero_slider", 
                             label = "Genes With at Least `X` Samples That Are Non-Zero:", 
                             min = 1, 
                             max = 69, 
                             value = 5
                 ),
                 actionButton(inputId = "reset_counts_main", 
                              label = "Reset Sliders")
               ),
               mainPanel( # Main panel for counts main-tab
                 tabsetPanel( # Sub-tabs panel for counts main-tab
                   tabPanel("Summary", # Summary sub-tab for counts main-tab
                            tags$h5("Please wait a few seconds for the summary to load."),
                            tableOutput("counts_summary_table")
                   ),
                   tabPanel("Diagnostic Plots", # Summary sub-tab for counts main-tab
                            tags$h5("Please wait a few seconds for the plots to load."),
                            tags$h5("Scatter plot of median count vs variance"),
                            plotOutput("plot_count_vs_variance"),
                            tags$h5("Scatter plot of median count vs number of zeros"),
                            plotOutput("plot_count_vs_zeros")
                   ),
                   # Tables sub-tab
                   tabPanel("Clustered Heatmap",
                            tags$div(
                              style = "margin-top: 5px;",
                              radioButtons(inputId = "count_heatmap_transformation", 
                                           label = "Select Whether to Log Transform Data:", 
                                           choices = c("TRUE", "FALSE"), 
                                           selected = "TRUE"
                              )
                            ),
                            actionButton(inputId = "count_heatmap_button", 
                                         label = "Display Heatmap"
                            ),
                            tags$h5("Please wait a few seconds for the plots to load."),
                            plotOutput("count_heatmap")
                   ),
                   # Plots sub-tab 
                   tabPanel("PCA", 
                            sidebarLayout( # Sidebar layout for counts main-tab
                              sidebarPanel(
                                radioButtons(inputId = "count_pca_x", 
                                             label = "Select which PC should be used for X-axis in PCA plot:", 
                                             choiceNames = list("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
                                             choiceValues = list("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
                                             selected = "PC1"
                                ),
                                radioButtons(inputId = "count_pca_y", 
                                             label = "Select which PC should be used for Y-axis in PCA plot:", 
                                             choiceNames = list("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
                                             choiceValues = list("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
                                             selected = "PC2"
                                ),
                                actionButton(inputId = "reset_count_pca_button", 
                                             label = "Default PCA")
                              ),
                              mainPanel(
                                plotOutput("count_pca_plot")
                              )
                            )
                   ) 
                 ) 
               ) 
             )
    ),
    
    # DE main-tab 
    tabPanel("DE", 
              tags$h3("Differential Expression Analysis"),
              tags$div(
                tags$h5("Visualize and search differential expression results")
              ),
             sidebarLayout( # Sidebar layout for DE main-tab
               sidebarPanel( # Sidebar panel for DE main-tab
                 fileInput(inputId = "de", # Input for uploading DE file 
                           label = "Please Upload Differential Expression File", 
                           multiple = FALSE, 
                           accept = ".csv",
                           buttonLabel = "Browse DEG File",
                           placeholder = "data/deg.csv"
                 ),
                 radioButtons("x_name", "Choose the column for the x-axis",
                              choices = c("baseMean" = "baseMean",
                                          "log2FoldChange" = "log2FoldChange",
                                          "lfcSE" = "lfcSE",
                                          "stat" = "stat",
                                          "pvalue" = "pvalue",
                                          "padj" = "padj"),
                              selected = "log2FoldChange"),
                 radioButtons("y_name", "Select Y Variable:",
                              choices = c("baseMean" = "baseMean",
                                          "log2FoldChange" = "log2FoldChange",
                                          "lfcSE" = "lfcSE",
                                          "stat" = "stat",
                                          "pvalue" = "pvalue",
                                          "padj" = "padj"),
                              selected = "padj"),
                 colourInput("base_color", "Base point color",
                             value = "#22577A"),
                 colourInput("highlight_color", "Highlight point color",
                             value = "#FFCF56"),
                 sliderInput("volcano_slider", "Select the magnitude of the p adjusted coloring:",
                             min = -300, max = 0, value = -10, step = 1),
                 actionButton(inputId = "reset_de_volcano", 
                              label = "Reset Slider"),
                 tags$hr(),
                 actionButton("plot_button", "Plot")
               ), 
               mainPanel( # Main panel for DE main-tab
                 tabsetPanel( # Sub-tabs panel for DE main-tab
                   tabPanel("Tables", # Tables sub-tab for DE main-tab
                            DTOutput("deg_table")
                   ), 
                   tabPanel("Volcano Plot", # Volcano plot sub-tab for DE main-tab
                            plotOutput("volcano_plot")
                   ),
                   tabPanel("Plot Table", # MA plot sub-tab for DE main-tab
                            DTOutput("volcano_plot_table")
                   ),
                 ) 
               ) 
             ) 
    ),
    
    # Individual Gene Expression main-tab 
    tabPanel("Individual Gene Expression", 
              tags$h3("Individual Gene Expression"), 
              tags$div(
                tags$h5("Visualize, search and plot individual gene expression")
              ),
             sidebarLayout(
               sidebarPanel(
                 fileInput(inputId = "sample_info",
                           label = "Upload Sample Information Matrix (CSV)",
                           multiple = FALSE,
                           accept = ".csv",
                           buttonLabel = "Browse .csv",
                           placeholder = "data/metadata.csv"
                 ),
                 fileInput(inputId = "norm_counts",
                           label = "Upload Normalized Counts Matrix (CSV)",
                           multiple = FALSE,
                           accept = ".csv",
                           buttonLabel = "Browse .csv",
                           placeholder = "data/norm_counts.csv"
                 ),
                 uiOutput("cat_field"),
                 uiOutput("gene_selector"),
                 radioButtons("plot_type", "Select plot type:",
                              choices = c("Bar plot" = "bar",
                                          "Box plot" = "box",
                                          "Violin plot" = "violin",
                                          "Beeswarm plot" = "beeswarm",
                                          "Scatter plot" = "scatter"),
                              selected = "bar"),
                 actionButton(inputId = "unique_plot_button", label = "Plot")
               ),
               mainPanel(
                 plotOutput("count_plot")
               )
             )
    )
  )
)


server <- function(input, output, session) {
  options(shiny.maxRequestSize=50*1024^2)
  # Reset sliders in counts main-tab
  # Existing Code
  observeEvent(input$reset_counts_main, {
    updateSliderInput(session = session, inputId = "count_percentile_slider", value = 80)
    updateSliderInput(session = session, inputId = "count_nonzero_slider", value = 5)
  })
  
  # 1. output functions for "Summary" sub-tab in "Samples" main-tab
  sample_data <- reactive({
    req(input$samples)
    metadata_raw <- rio::import(input$samples$datapath)
    # Generate summary paragraph based on sample metadata
    metadata_info_df <- metadata_raw %>%
      dplyr::select(geo_accession, publish_date = status, organ = source_name_ch1, species = organism_ch1,
                    tissue = `tissue:ch1`, seq = library_strategy, lib = library_source, 
                    lib_selection = library_selection, instrument = instrument_model, molecule = molecule_ch1, 
                    extract_prot1 = extract_protocol_ch1, extract_prot2 = extract_protocol_ch1.1, 
                    tax_id = taxid_ch1, data_proc1 = data_processing, data_proc2 = data_processing.1,       
                    data_proc3 = data_processing.2, data_proc4 = data_processing.3, data_proc5 = data_processing.4,      
                    data_proc6 = data_processing.5, data_proc7 = data_processing.6, data_proc8 = data_processing.7,       
                    platform_id, contact_name, contact_email, contact_dept = contact_department,      
                    contact_inst = contact_institute, contact_address, contact_city, contact_state,           
                    contact_zip = `contact_zip/postal_code`, contact_country)
    metadata_info <- metadata_info_df %>%
      dplyr::distinct() %>%
      tibble::remove_rownames() %>%
      dplyr::mutate(metadata_info = paste0(
        "<b>Data Publication:</b> The data was published on ", publish_date, 
        " with ", geo_accession, 
        " and involved sequencing the ", molecule,
        " of ", organ,
        " from ", species,
        " using the ", instrument,
        " platform.<br><br>",
        "<b>Tissue and Library:</b> The tissue type was ", tissue,
        " and the library source was ", lib,
        " with the library selection method of ", lib_selection,
        ". The extraction protocols were ", extract_prot1,
        " and ", extract_prot2,
        ". The taxonomic ID was ", tax_id,
        ".<br><br>",
        "<b>Data Processing:</b> The data processing steps were as follows:<br>",
        "- ", data_proc1,
        "<br>- ", data_proc2,
        "<br>- ", data_proc3,
        "<br>- ", data_proc4,
        "<br>- ", data_proc5,
        "<br>- ", data_proc6,
        "<br>- ", data_proc7,
        "<br>- ", paste("Genome build:", data_proc8),
        ".<br><br>",
        "<b>Platform ID:</b> The platform ID was ", platform_id,
        ".<br><br>",
        "<b>Contact Information:</b> The contact information for this project was:<br>",
        "- Name: ", contact_name,
        "<br>- Email: ", contact_email,
        "<br>- Department: ", contact_dept,
        "<br>- Institute: ", contact_inst,
        "<br>- Address: ", paste(contact_address, ",", contact_city, ",", contact_state, ",", contact_zip, ",", contact_country),
        ".")) %>%
      head(1) %>%
      dplyr::pull(metadata_info)
    # Acquire the summary table of sample metadata
    metadata <- metadata_raw %>%
      dplyr::select(Diagnosis = `diagnosis:ch1`, 
                    Organ = source_name_ch1, 
                    Species = organism_ch1,
                    Tissue = `tissue:ch1`, 
                    RIN = `rin:ch1`, 
                    PMI = `pmi:ch1`,
                    Sequence_Reads = `mrna-seq reads:ch1`,
                    Age_of_Death = `age of death:ch1`,       
                    Age_of_Onset = `age of onset:ch1`,        
                    CAG = `cag:ch1`,                
                    Duration = `duration:ch1`,           
                    HV_Cortical_Score = `h-v cortical score:ch1`,  
                    HV_Striatial_Score = `h-v striatal score:ch1`, 
                    Vonsattel_Grade = `vonsattel grade:ch1`,
                    Channel_Count = `channel_count`) %>%
      dplyr::mutate(across(-c(Diagnosis, Organ, Species, Tissue), as.numeric))
    metadata_summary <- metadata %>%
      dplyr::reframe(Type = map(., typeof),
                     Mean = map(., function(x) ifelse(is.numeric(x), mean(na.omit(x)), NA)),
                     SD = map(., function(x) ifelse(is.numeric(x), sd(na.omit(x)), NA)),
                     Distinct_Values = map(., function(x) ifelse(is.numeric(x), NA, paste(unique(na.omit(x)), collapse = ", ")))) %>%
      dplyr::mutate(`Column Name` = colnames(metadata)) %>%
      dplyr::relocate(`Column Name`, .before = Type)
    metadata <- dplyr::mutate(metadata, geo_accession = metadata_raw$title)
    plot_data <- dplyr::select(metadata, 
                               Samples = geo_accession, PMI, Diagnosis, `RIN Score` = RIN, 
                               `mRNA Reads` = Sequence_Reads, `Age of Death` = Age_of_Death)
    return(list(summary_table = metadata_summary, summary_paragraph = metadata_info, 
                metadata = metadata, metadata_info_df = metadata_info_df,
                metadata_raw = metadata_raw, plot_data = plot_data))
  })
  # Rendering summary table and paragraph
  output$sample_summary <- renderTable({
    sample_data()$summary_table
  })
  output$sample_info <- renderText({
    sample_data()$summary_paragraph
  })
  
  
  # 2. output functions for "Tables" sub-tab in "Samples" main-tab (both short and long tables)
  # Reactive value to track the state of the table to be displayed
  display_table <- reactiveVal("sample_short_table")
  observeEvent(input$sample_simple_table_display_button, {
    display_table("sample_short_table")
  })
  observeEvent(input$sample_table_display_button, {
    display_table("sample_table")
  })
  # Displaying short tables
  output$sample_short_table <- renderDT({
    metadata_raw <- rio::import(input$samples$datapath)
    table <- dplyr::select(metadata_raw, c("geo_accession", "title", input$sample_display_col))
    # Rename columns
    old_names <- c("diagnosis:ch1", "rin:ch1", "pmi:ch1", "mrna-seq reads:ch1", "age of death:ch1",
                   "age of onset:ch1", "cag:ch1", "duration:ch1", "h-v cortical score:ch1", "h-v striatal score:ch1", "vonsattel grade:ch1")
    new_names <- c("Diagnosis", "RIN", "PMI", "Sequence_Reads", "Age_of_Death", "Age_of_Onset", "CAG", "Duration", "HV_Cortical_Score",
                   "HV_Striatial_Score", "Vonsattel_Grade")
    new_names <- new_names[old_names %in% input$sample_display_col]
    colnames(table) <- c("GEO", "Samples", new_names)
    table <- dplyr::arrange(table, desc(Samples))
  })
  # Displaying long tables
  output$sample_table <- renderDT({
    req(input$sample_table_display_button)
    metadata <- sample_data()$metadata
    metadata_info_df <- sample_data()$metadata_info_df
    table <- dplyr::right_join(metadata_info_df, metadata, by = "geo_accession")
  })
  output$sample_table_display <- renderUI({
    if (display_table() == "sample_short_table") {
      DTOutput("sample_short_table")
    } else {
      DTOutput("sample_table")
    }
  })
  
  
  # 3. output plot functions for "Plots" sub-tab in "Samples" main-tab
  # Histogram and density plot
  display_plot <- reactiveVal("histogram")
  observeEvent(input$plot_histogram, {
    display_plot("histogram")
  })
  observeEvent(input$plot_density, {
    display_plot("density")
  })
  # Render Plots
  output$histogram_plot <- renderPlot({
    req(input$sample_plot_dropdown)
    var_name <- input$sample_plot_dropdown
    data <- sample_data()$plot_data
    ggplot(data, aes(x = !!sym(var_name))) +
      geom_histogram() +
      theme_classic()
  })
  output$density_plot <- renderPlot({
    req(input$sample_plot_dropdown)
    var_name <- input$sample_plot_dropdown
    data <- sample_data()$plot_data
    ggplot(data, aes(x = !!sym(var_name))) +
      geom_density() +
      theme_classic()
  })
  output$sample_plot_display <- renderUI({
    if (display_plot() == "histogram") {
      plotOutput("histogram_plot")
    } else if (display_plot() == "density") {
      plotOutput("density_plot")
    }
  })
  
  
  # 4. plot functions for "Summary" sub-tab in "Counts" main-tab
  # Filter the data based on the sliders
  filtered_counts <- reactive({
    req(input$counts)
    counts_data <- rio::import(input$counts$datapath)
    
    # Calculate the variance for each gene (row)
    variances <- apply(counts_data[, -1], 1, var)
    # Calculate the quantile based on the selected percentile
    variance_threshold <- quantile(variances, input$count_percentile_slider / 100)
    # Filter genes based on the variance threshold
    filtered_by_variance <- counts_data[variances >= variance_threshold, ]
    # Count the number of non-zero samples for each gene (row)
    nonzero_counts <- rowSums(filtered_by_variance != 0)
    # Filter genes based on the non-zero samples threshold
    filtered_by_nonzero <- filtered_by_variance[nonzero_counts >= input$count_nonzero_slider, ]
    filtered_data <- filtered_by_nonzero
    
    num_samples <- ncol(counts_data) - 1
    num_genes <- nrow(counts_data)
    num_passed_filter <- nrow(filtered_data)
    percent_passed_filter <- (num_passed_filter / num_genes) * 100
    num_not_passed_filter <- num_genes - num_passed_filter
    percent_not_passed_filter <- (num_not_passed_filter / num_genes) * 100
    
    summary_table <-   tibble::tribble(
      ~Measure, ~Summary,
      "Number of Samples", num_samples,
      "Number of Genes", num_genes,
      "Number of Genes Passed Filter", num_passed_filter,
      "% Passed Filter", percent_passed_filter,
      "Number of Genes Not Passed Filter", num_not_passed_filter,
      "% Not Passed Filter", percent_not_passed_filter
    )
    return(list(counts_data = counts_data, filtered_by_nonzero = filtered_by_nonzero, 
                summary_table = summary_table, variances = variances, variance_threshold = variance_threshold))
  })
  
  # Calculate and output the summary table
  output$counts_summary_table <- renderTable({
    filtered_counts()$summary_table
  })
  
  
  # 5. plot functions for "Diagnostic Plots" sub-tab in "Counts" main-tab
  # Diagnostic plot: Median count vs. variance
  output$plot_count_vs_variance <- renderPlot({
    counts_data <- filtered_counts()$counts_data
    variances <- filtered_counts()$variances
    variance_threshold <- filtered_counts()$variance_threshold
    
    medians <- apply(counts_data[, -1], 1, median)
    df <- data.frame(medians, variances) %>%
      dplyr::mutate(Threshold = ifelse(variances < variance_threshold, "Failed Threshold", "Passed Threshold"))
    
    ggplot(df, aes(x = medians, y = variances, color = Threshold)) +
      geom_point() +
      scale_color_manual(values = c("#90EE90", "black")) +
      scale_x_log10() +
      scale_y_log10() +
      labs(x = "Log 10 Median Count", y = "Log 10 Variance") +
      theme_classic() +
      theme(legend.position = "bottom")
  })
  
  # Diagnostic plot: Median count vs. number of zeros
  output$plot_count_vs_zeros <- renderPlot({
    counts_data <- filtered_counts()$counts_data
    
    medians <- apply(counts_data[, -1], 1, median)
    num_zeros <- apply(counts_data[, -1], 1, function(x) sum(x == 0))
    df <- data.frame(medians, num_zeros) %>%
      dplyr::mutate(Threshold = ifelse(num_zeros > input$count_nonzero_slider, "Failed Filter", "Passed Filter"))
    
    # Create a scatter plot of median count vs. number of zeros
    plot <- ggplot(df, aes(x = medians, y = num_zeros, color = Threshold)) +
      geom_point() +
      scale_color_manual(values = c("#90EE90", "black")) +
      scale_x_log10() +
      scale_y_log10() +
      ylim(0, 70) +
      labs(x = "Log 10 Median Count", y = "Log 10 Number of Samples with Zero Counts") +
      theme_classic() +
      theme(legend.position = "bottom")
    plot
  })
  
  # 6. plot function for "Clustered Heatmap" sub-tab in "Counts" main-tab
  count_heatmap_radio <- reactive({
    if (is.null(input$count_heatmap_transformation)) {
      return("TRUE")
    }
    input$count_heatmap_transformation
  })
  output$count_heatmap <- renderPlot({
    req(input$count_heatmap_button)
    data <- filtered_counts()$filtered_by_nonzero %>%
      data.frame() %>%
      tibble::remove_rownames() %>%
      dplyr::select(-gene) %>%
      as.matrix()
    if(count_heatmap_radio() == "TRUE") {
      data <- log(data + 1)
    } 
    isolate({
      heatmap.2(data)
    })
  })
  
  # 7. plot function for PCA in "Counts" main-tab
  observeEvent(input$reset_count_pca_button, {
    updateRadioButtons(session = session, inputId = "count_pca_x", selected = "PC1")
    updateRadioButtons(session = session, inputId = "count_pca_y", selected = "PC2")
  })
  
  output$count_pca_plot <- renderPlot({
    data <- filtered_counts()$filtered_by_nonzero %>%
      data.frame() %>%
      tibble::remove_rownames() %>%
      dplyr::select(-gene) %>%
      as.matrix()
    pca <- prcomp(data)
    pca_df <- as.data.frame(pca$x)
    ggplot(pca_df, aes(x=!!sym(input$count_pca_x), y=!!sym(input$count_pca_y))) +
      geom_point() +
      xlab(input$count_pca_x) +
      ylab(input$count_pca_y) +
      ggtitle("PCA Plot") +
      theme_classic()
  })
  
  
  # 8. a plot/table function for a DEG table in the "DE" main-tab
  observeEvent(input$reset_de_volcano, {
    updateSliderInput(session, "volcano_slider", value = -10)
  })
  
  de_data <- reactive({
    req(input$de)
    deg <- rio::import(input$de$datapath)
    return(deg)
  })
  
  # Display the DE data table
  output$deg_table <- renderDT({
    deg <- de_data()
  })
  
  # Create the volcano plot
  output$volcano_plot <- renderPlot({
    req(input$plot_button)
    deg <- de_data()
    isolate({
      df <- de_data()
      ggplot(df, aes(x = !!sym(input$x_name), y = -log10(!!sym(input$y_name)))) +                
        geom_point(aes(color = (df[[input$y_name]] <= 10^input$volcano_slider)), size = 1) +
        scale_color_manual(values = c(input$base_color, input$highlight_color)) +
        labs(x = input$x_name, y = input$y_name) + 
        theme_classic() + 
        theme(legend.position = "bottom")
    })
  })
  
  # Display the DE data table with colors
  output$volcano_plot_table <- renderDT({
    req(input$plot_button)
    df <- de_data() %>%
      dplyr::filter(!!sym(input$y_name) <= 10^input$volcano_slider)
  })
  
  
  # 9. Individual Gene Expression:
  norm_counts <- reactive({
    if(is.null(input$norm_counts)){
      counts <- filtered_counts()$counts_data
    } else{
      counts <- rio::import(input$norm_counts$datapath)
    }
    counts <- counts
  })
  
  sample_info <- reactive({
    if(is.null(input$sample_info)){
      raw_metadata <- sample_data()$metadata_raw
    } else{
      raw_metadata <- rio::import(input$sample_info$datapath)
    } 
    metadata <- raw_metadata %>%
      dplyr::select(Diagnosis = `diagnosis:ch1`,
                    Samples = title,
                    RIN = `rin:ch1`, 
                    PMI = `pmi:ch1`,
                    Sequence_Reads = `mrna-seq reads:ch1`,
                    Age_of_Death = `age of death:ch1`) %>%
      dplyr::mutate(across(-c(Diagnosis, Samples), as.numeric),
                    Diagnosis = as.factor(Diagnosis))
    metadata
  })
  
  # Update categorical field
  output$cat_field <- renderUI({
    si <- sample_info()
    cat_fields <- colnames(si)
    selectInput("selected_cat_field", "Choose categorical field:", choices = cat_fields)
  })
  output$gene_selector <- renderUI({
    genes <- norm_counts()$gene
    selectInput("selected_gene", "Select a gene:", choices = genes)
  })
  
  
  output$count_plot <- renderPlot({
    req(input$unique_plot_button)
    req(norm_counts(), sample_info())
    isolate({
      counts <- norm_counts()
      metadata <- sample_info()
      # Merge the data based on input
      merged_data <- counts %>%
        dplyr::filter(gene == input$selected_gene) %>% 
        dplyr::select(-gene) %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("Samples") %>% 
        dplyr::left_join(metadata, by = "Samples") %>% 
        dplyr::rename(count = V1)
      
      
      # Create the plot based on the selected plot type
      switch(input$plot_type,
             bar = ggplot(merged_data, aes_string(x = input$selected_cat_field, y = "count")) +
               geom_bar(stat = "identity", position = "dodge") +
               labs(title = "Bar Plot", x = input$selected_cat_field, y = "Normalized Count"),
             box = ggplot(merged_data, aes_string(x = input$selected_cat_field, y = "count")) +
               geom_boxplot() +
               labs(title = "Box Plot", x = input$selected_cat_field, y = "Normalized Count"),
             violin = ggplot(merged_data, aes_string(x = input$selected_cat_field, y = "count")) +
               geom_violin() +
               labs(title = "Violin Plot", x = input$selected_cat_field, y = "Normalized Count"),
             beeswarm = ggplot(merged_data, aes_string(x = input$selected_cat_field, y = "count")) +
               geom_beeswarm(width = 0.4) +
               labs(title = "Beeswarm Plot", x = input$selected_cat_field, y = "Normalized Count"),
             scatter = ggplot(merged_data, aes_string(x = input$selected_cat_field, y = "count")) +
               geom_jitter() +
               labs(title = "Scatter Plot", x = input$selected_cat_field, y = "Normalized Count")
      ) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
  })
  
  
  
  output$testTest <- renderPrint({
    test1 <- input$count_percentile_slider
    test2 <- input$count_nonzero_slider
    paste(test1, test2)
  })
}

shinyApp(ui, server)