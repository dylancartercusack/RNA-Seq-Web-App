
# Load packages ----------------------------------------------------------------------------

library(shiny)
library(ggplot2)
library(tidyverse)
library(corrr)
library(shinythemes)
library(plotly)
library(thematic)
#library(EnhancedVolcano)
library(forcats)
library(ggrepel)
library(dplyr)
library(ggdark)
library(rsconnect)
library(ggtext)

# Load and prepare data --------------------------------------------------------------------
# 
print("Starting up....")
# setwd("/Volumes/tri_services/Facilities/Group_Data/MR_Groups/RG24. Macrophage Biology - DH/Dylan Carter-Cusack/Bioinformatics/rat_tissue_seq_v3/web_app/final")
# load("./web_app_data_final.RData")



# Data preparation - don't need this every time -------------------------------------------
# Read in expression file
# raw_data <- read.csv("/Volumes/tri_services/Facilities/Group_Data/MR_Groups/RG24. Macrophage Biology - DH/Dylan Carter-Cusack/Bioinformatics/rat_tissue_seq_v3/web_app/test2/master_tpm_collated_12112023.csv")
# # #
# # #  # Separate out metadata on its own
# metadata <- filter(raw_data, row_number() %in% c(1:5)) %>%
#   subset(select = -c(gene_name, description)) %>%
#   column_to_rownames("unique_gene_id") %>%
# t() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "sample")
# # #
# # #  # Separate out expression data on its own
# cleaned_expression <- filter(raw_data, !row_number() %in% c(1:5)) %>%
#   subset(select = -c(gene_name, description)) %>%
#   pivot_longer(cols = !1, names_to = "sample", values_to = "tpm")
# 
# # #  # Join expression data and metadata after they have been converted to long format
# final_data <- cleaned_expression %>%
#   inner_join(metadata, by = "sample") %>%
#   transform(tpm = as.numeric(tpm)) %>%
#   transform(sample = as.factor(sample))
# 
# # #  # Calculate mean and standard deviation for each gene by tissue type and condition
# final_data <- final_data %>% group_by(unique_gene_id, tissue, condition) %>%
#   mutate(mean_tpm = mean(tpm), sd_tpm = sd(tpm)) %>%
#   ungroup()
# #
# gene_list <- filter(raw_data, !row_number() %in% c(1:4)) %>%
#   subset(select = 1)
# # # 
# full_DE_list <- read.csv("/Volumes/tri_services/Facilities/Group_Data/MR_Groups/RG24. Macrophage Biology - DH/Dylan Carter-Cusack/Bioinformatics/rat_tissue_seq_v3/web_app/test2/full_DE_list.csv")

#save.image("/Volumes/tri_services/Facilities/Group_Data/MR_Groups/RG24. Macrophage Biology - DH/Dylan Carter-Cusack/Bioinformatics/rat_tissue_seq_v3/web_app/final/web_app_data_final.RData")


# # Create full list of DE genes from each tissue
# wat <- read.csv("./data/DE_lists/wat.csv")
# full_DE_list <- bind_rows(lst(adrenal, bone, cerebellum, colon, cortex, diaphragm, dLN, duodenum, heart, hippocampus, hypothalamus, ileum, kidney, liver, lung, mLN, olfactorybulb, ovary, pituitary, skmuscle, testes, thymus, tibialisant, wat), .id = "id") %>%
#   dplyr::rename(tissue = id)


# Define UI ---------------------------------------------------------------------------------

ui <- fluidPage(
  #shinythemes::themeSelector(),
  theme = shinytheme("slate"),
  navbarPage("Rat Atlas - 3 week WT and Csf1rko Tissue Bulk RNA-Seq", windowTitle = "Rat Atlas v2"),
  
  fluidRow(
    column(2,
      # Tissue selection
      selectizeInput(
        inputId = "tissue_selection",
        label = "Select a tissue:",
        choices = list("All Tissues" = "all_tissues","Adrenal Gland" = "adrenal", "Colon" = "colon", "Ileum" = "ileum", "Duodenum" = "duodenum", "Heart" = "heart", "Diaphragm" = "diaphragm", "Tibialis Anterior" = "tibialisant", "Quadriceps" = "quadriceps", "Thymus" = "thymus", "Lymph Nodes" = "LN", "Kidney" = "kidney", "Lung" = "lung", "Ovary" = "ovary", "Testes" = "testes", "White Adipose Tissue" = "wat", "Liver" = "liver", "Cerebellum" = "cerebellum", "Cortex" = "cortex", "Hippocampus" = "hippocampus", "Hypothalamus" = "hypothalamus", "Olfactory Bulb" = "olfactorybulb", "Pituitary" = "pituitary", "Transplanted Lung" = "lttx_lung", "Transplanted Liver" = "lttx_liver", "Transplanted Peritoneal Cavity" = "lttx_PT", "Transplanted Brain" = "lttx_brain"),
        selected = NULL,
        multiple = FALSE
      ),
      
    ),
    # Text input for gene and submit button to inititate search
    column(2,
      selectizeInput(
        inputId = "unique_gene_id",
        label = "Gene Name:",
        choices = NULL,
        #options = list(create = TRUE),
        selected = NULL
    )),
    # Submit Button
    column(2,
    actionButton(
      inputId = "submit_button",
      label = "Submit",
      class = "btn-primary"# Align button with the input box
    )), 
    
    # Add the note at the top of the page
    column(2),
    column(4,
    tags$div(
      style = "background-color: #aba8a6; padding: 10px; border: 1px solid #ddd; margin-bottom: 10px;",
      tags$p(
        style = "font-size: 14px; color: #555;",
        "Note: These plots are interactive. Click and drag anywhere on the plot to zoom in. To change only the y-axis, click and drag straight up or down. To reset the zoom, either double click anywhere on the plot, or use the icons in the top right of the plot."
      )
    )),
    
    ),
  #Output: Show scatterplot - overflow-y: scroll; position: relative"
  fluidRow(
    column(12,
         
         tabsetPanel( 
           tabPanel(title = "Barplot",
                    (div(style = "max-height:300px;",
                         uiOutput("barplot.ui")
                    ))),
           
           tabPanel(title = "Volcano Plot",
                    fluidRow(column(5, div(plotlyOutput(
                      outputId = "volcano_plot",
                      height = "800px"
                    ))),
                    column(7, div(dataTableOutput('DElist_subset')
                    ))
                    ),
  )
)

)
)
)

# Define server -----------------------------------------------------------------------------

server <- function(input, output, session) {
  #thematic::thematic_shiny()

  
  ##### --------------
  # Add functionality to go to the Ensembl page for the selected gene
  # EIther use the API or just find how the web link changes for each gene and then inject the user input gene into that
  ##### --------------
  
  
  # Use server side selectize to create gene selection list
  updateSelectizeInput(session = getDefaultReactiveDomain(),
                       "unique_gene_id",
                       choices = gene_list$unique_gene_id,
                       selected = "Csf1r",
                       server = TRUE)
  
  
  # Subset dataset to only include required genes and samples
  df_subset <- reactive({
    if (input$tissue_selection == "all_tissues") { # If the tissue selection is all tissues, then just filter down to the selected gene
      filter(final_data, unique_gene_id == input$unique_gene_id)
    } else { # If a specific tissue is selected, then filter down to just that tissue, and then filter down to selected gene
      filter(final_data, grepl(input$tissue_selection, sample)) %>%
        filter(unique_gene_id == input$unique_gene_id)
    }
  })
  
    plot_height <- 900
    # plot_height <- bindEvent(
    # reactive({500+30*nrow(df_subset())}),
    # input$submit_button)
  
    
    group_order <- c("adrenal_ko", "adrenal_wt", "colon_ko", "colon_wt", "ileum_ko", "ileum_wt", "duodenum_ko", "duodenum_wt", "kidney_ko", "kidney_wt", "lung_ko", "lung_wt", "ovary_ko", "ovary_wt", "teste_ko", "teste_wt", "liver_ko", "liver_wt", "bone_ko", "bone_wt", "wat_ko", "wat_wt", "thymus_ko", "thymus_wt", "LN_ko", "LN_wt", "heart_ko", "heart_wt", "diaphragm_ko", "diaphragm_wt", "tibialisant_ko", "tibialisant_wt", "quadriceps_ko", "quadriceps_wt", "cerebellum_ko", "cerebellum_wt", "cortex_ko", "cortex_wt", "hippocampus_ko", "hippocampus_wt", "hypothalamus_ko", "hypothalamus_wt", "olfactory_bulb_ko", "olfactory_bulb_wt", "pituitary_ko", "pituitary_wt", "striatum_ko", "striatum_wt", "lttx_brain_ko", "lttx_brain_wt","lttx_lung_ko", "lttx_lung_wt", "lttx_liver_ko", "lttx_liver_wt", "lttx_PT_ko", "lttx_PT_wt")
    graph_labels <- ifelse(final_data$condition == "ko", "red", "blue")

  # Barplot of selected gene
  output$barplot <- bindEvent(
    renderPlotly({
      #plot_height <- 700+100*nrow(df_subset())
      plot <-
        ggplot(df_subset(), aes(x = factor(tissue_condition, level = group_order), y = tpm, colour = tissue, group = tissue_condition)) +
        geom_point(aes(colour = tissue), position = position_jitterdodge(0.4)) +
        geom_boxplot(aes(colour = tissue), alpha = 0.2) +
        dark_theme_gray() +
        #coord_flip() +
        theme(axis.text.x = element_text(angle=90, size = 15), axis.title.x = element_blank(), axis.text.y = element_text(size = 15),
              legend.position = "none", legend.direction = "horizontal", text = element_text(family = "Helvetica", size = 15))
        
       plot %>%
        ggplotly(tooltip = "all", dynamicTicks = "y") %>% 
         # rangeslider() %>% 
        layout(legend = list(orientation = "h", x = 0.45, y = 1))
  }),
  input$submit_button)

  output$barplot.ui <- renderUI({
    plotlyOutput("barplot", height = plot_height) #height = plot_height())
  })


  output$volcano_plot <- 
    bindEvent(
    renderPlotly({
        # Take input from sliders
        volcano_logFC_cutoff <- 1 #input$logFC_slider
        volcano_padj_cutoff <- 0.05 #input$padj_slider
        
        # Subset full_DE_list to only include the selected tissue
        DElist_subset <- filter(full_DE_list, tissue == input$tissue_selection) %>% 
          select(!c(tissue, baseMean, lfcSE, stat, gene_id))
        # Create the DE table to output
        output$DElist_subset <- renderDataTable(DElist_subset, options = list(pageLength = 15, autoWidth = TRUE))
        
        
        DElist_subset_plot <- DElist_subset
        DElist_subset_plot$diffexpressed <- "NO"
        DElist_subset_plot$diffexpressed[DElist_subset_plot$log2FoldChange > volcano_logFC_cutoff & DElist_subset_plot$padj < volcano_padj_cutoff] <- "UP"
        DElist_subset_plot$diffexpressed[DElist_subset_plot$log2FoldChange < -(volcano_logFC_cutoff) & DElist_subset_plot$padj < volcano_padj_cutoff] <- "DOWN"
        DElist_subset_plot$delabel <- NA
        DElist_subset_plot$delabel[DElist_subset_plot$diffexpressed != "NO"] <- DElist_subset_plot$unique_gene_id[DElist_subset_plot$diffexpressed != "NO"]
        DElist_subset_plot$gene_of_interest <- "NO"
        DElist_subset_plot$gene_of_interest[DElist_subset_plot$unique_gene_id == input$unique_gene_id] <- "HIGHLIGHT"     #DElist_subset_plot$unique_gene_id[DElist_subset_plot$unique_gene_id == input$unique_gene_id]
        

      # Selected
      # volcano_colours <- c("blue", "red", "black")
      # names(volcano_colours) <- c("DOWN", "UP", "NO")
      
      
      # Draw volcano plot
      volcano_plot <- ggplot(DElist_subset_plot, 
                             aes(x = log2FoldChange, y = -log10(padj), label = unique_gene_id,
                                 colour = interaction(diffexpressed, gene_of_interest, sep = "_"))) +
        geom_point() +
        geom_vline(xintercept = c(-(volcano_logFC_cutoff), volcano_logFC_cutoff), col = "red") +
        geom_hline(yintercept = -log10(volcano_padj_cutoff), col = "red") +
        dark_theme_gray() +
        theme(legend.position = "none", text = element_text(family = "Helvetica", size = 15))
        
      
      #+
        # scale_colour_manual(values = volcano_colours)
      
      
      
      
      volcano_plot %>%
        ggplotly(tooltip = c("unique_gene_id", "log2FoldChange", "padj", "-log10(padj)"))# %>%
        # add_trace(
        #   type = "scatter",
        #   hovertemplate = paste('%{unique_gene_id')
      # )
        

    }),
    input$submit_button
    )

  
  
}


# Create a Shiny app object ----------------------------------------------------
shinyApp(ui = ui, server = server)








# volcano_padj_cutoff <- 10e-5
# volcano_logFC_cutoff <- 1.5
# DElist_subset()$diffexpressed <- "NO"
# DElist_subset()$diffexpressed[DElist_subset()$log2FoldChange > volcano_logFC_cutoff & DElist_subset()$padj < volcano_padj_cutoff] <- "UP"
# DElist_subset()$diffexpressed[DElist_subset()$log2FoldChange < -(volcano_logFC_cutoff) & DElist_subset()$padj < volcano_padj_cutoff] <- "DOWN"
# DElist_subset()$delabel <- NA
# DElist_subset()$delabel[DElist_subset()$diffexpressed != "NO"] <- DElist_subset()$unique_gene_id[DElist_subset()$diffexpressed != "NO"]
# volcano_colours <- c("blue", "red", "black")
# names(volcano_colours) <- c("DOWN", "UP", "NO")
# 
# volcano_plot <- ggplot(DElist_subset(), aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
#   geom_point() +
#   geom_text_repel() +
#   geom_vline(xintercept = c(-(volcano_logFC_cutoff), volcano_logFC_cutoff), col = "red") +
#   geom_hline(yintercept = -log10(volcano_padj_cutoff), col = "red") +
#   scale_colour_manual(values = volcano_colours)
# volcano_plot
# 
# 
# 
# 
# 
# 
#








# 
# output$volcano_plot <- 
#   bindEvent(
#     renderPlotly({
#       
#       volcano_logFC_cutoff <- input$logFC_slider
#       volcano_padj_cutoff <- input$padj_slider
#       
#       
#       DElist_subset()$diffexpressed <- "NO"
#       DElist_subset()$diffexpressed[DElist_subset()$log2FoldChange > volcano_logFC_cutoff & DElist_subset()$padj < volcano_padj_cutoff] <- "UP"
#       DElist_subset()$diffexpressed[DElist_subset()$log2FoldChange < -(volcano_logFC_cutoff) & DElist_subset()$padj < volcano_padj_cutoff] <- "DOWN"
#       DElist_subset()$delabel <- NA
#       DElist_subset()$delabel[DElist_subset()$diffexpressed != "NO"] <- DElist_subset()$unique_gene_id[DElist_subset()$diffexpressed != "NO"]
#       DElist_subset()$gene_of_interest <- NA
#       reactive({DElist_subset()$gene_of_interest[DElist_subset()$unique_gene_id == input$unique_gene_id] <- DElist_subset()$unique_gene_id[DElist_subset()$unique_gene_id == input$unique_gene_id]
#       })
#       volcano_colours <- c("blue", "red", "black")
#       names(volcano_colours) <- c("DOWN", "UP", "NO")





