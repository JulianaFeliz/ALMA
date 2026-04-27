#!/usr/bin/env Rscript

# =============================================================================
# TraitAM Interactive Dashboard Launcher
# =============================================================================
# Description: Launch Shiny dashboard for exploring TraitAM results
# Author: Seqera AI
# Date: 2026-04-17
# Version: 1.0.0
# =============================================================================

suppressPackageStartupMessages({
    library(shiny)
    library(shinydashboard)
    library(tidyverse)
    library(DT)
    library(plotly)
    library(ggpubr)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

results_dir <- ifelse(length(args) >= 1, args[1], "results")

cat("\n=============================================================================\n")
cat("TraitAM Interactive Dashboard\n")
cat("=============================================================================\n\n")
cat("Loading results from:", results_dir, "\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("📂 Loading data files...\n")

# Find all sample results
sample_dirs <- list.dirs(file.path(results_dir, "02_trait_data"), 
                         recursive = FALSE, full.names = FALSE)

# If no subdirectories, look for files directly
trait_files <- list.files(
    file.path(results_dir, "02_trait_data"),
    pattern = "_traits_abundance\\.tsv$",
    full.names = TRUE
)

if (length(trait_files) == 0) {
    stop("No trait data files found in ", file.path(results_dir, "02_trait_data"))
}

# Extract sample IDs from filenames
sample_ids <- basename(trait_files) %>%
    str_remove("_traits_abundance\\.tsv$")

cat("  Found", length(sample_ids), "sample(s):", paste(sample_ids, collapse = ", "), "\n\n")

# Load all data
all_traits <- list()
all_fd <- list()
all_cwm <- list()

for (sample_id in sample_ids) {
    # Traits
    trait_file <- file.path(results_dir, "02_trait_data", 
                            paste0(sample_id, "_traits_abundance.tsv"))
    if (file.exists(trait_file)) {
        all_traits[[sample_id]] <- read_tsv(trait_file, show_col_types = FALSE) %>%
            mutate(sample_id = sample_id)
    }
    
    # FD metrics
    fd_file <- file.path(results_dir, "03_functional_diversity", 
                         paste0(sample_id, "_FD_metrics.tsv"))
    if (file.exists(fd_file)) {
        all_fd[[sample_id]] <- read_tsv(fd_file, show_col_types = FALSE)
    }
    
    # CWM
    cwm_file <- file.path(results_dir, "03_functional_diversity", 
                          paste0(sample_id, "_CWM.tsv"))
    if (file.exists(cwm_file)) {
        all_cwm[[sample_id]] <- read_tsv(cwm_file, show_col_types = FALSE)
    }
}

# Combine data
traits_combined <- bind_rows(all_traits)
fd_combined <- bind_rows(all_fd)
cwm_combined <- bind_rows(all_cwm)

cat("✓ Data loaded successfully\n")
cat("  Total records:", nrow(traits_combined), "\n")
cat("  Unique species:", n_distinct(traits_combined$species_name), "\n\n")

# =============================================================================
# SHINY APP UI
# =============================================================================

ui <- dashboardPage(
    skin = "blue",
    
    ## HEADER
    dashboardHeader(
        title = "TraitAM Analysis Dashboard",
        titleWidth = 300
    ),
    
    ## SIDEBAR
    dashboardSidebar(
        width = 300,
        sidebarMenu(
            menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
            menuItem("Species Explorer", tabName = "species", icon = icon("dna")),
            menuItem("Trait Analysis", tabName = "traits", icon = icon("chart-bar")),
            menuItem("Functional Diversity", tabName = "fd", icon = icon("layer-group")),
            menuItem("Community Comparison", tabName = "compare", icon = icon("balance-scale")),
            menuItem("Data Tables", tabName = "tables", icon = icon("table"))
        ),
        
        hr(),
        
        # Sample selector (if multiple samples)
        if (length(sample_ids) > 1) {
            selectInput(
                "sample_select",
                "Select Sample:",
                choices = sample_ids,
                selected = sample_ids[1]
            )
        } else {
            p(style = "padding: 10px;", strong("Sample:"), sample_ids[1])
        },
        
        hr(),
        
        # Trait selector for plots
        selectInput(
            "trait_select",
            "Select Trait:",
            choices = c(
                "Spore Volume" = "Spore_vol_um3",
                "Aspect Ratio" = "Aspect_ratio",
                "Color Score" = "Color_score",
                "Wall Investment" = "Wall_investment",
                "Ornamentation Height" = "Ornamentation_height_um"
            ),
            selected = "Spore_vol_um3"
        ),
        
        hr(),
        
        p(style = "padding: 10px; font-size: 11px;",
          "TraitAM Dashboard v1.0.0",
          br(),
          "Seqera AI © 2026")
    ),
    
    ## BODY
    dashboardBody(
        tabItems(
            # TAB 1: Overview
            tabItem(
                tabName = "overview",
                h2("Analysis Overview"),
                
                fluidRow(
                    valueBoxOutput("box_n_species"),
                    valueBoxOutput("box_n_families"),
                    valueBoxOutput("box_total_abundance")
                ),
                
                fluidRow(
                    box(
                        title = "Functional Diversity Metrics",
                        width = 6,
                        status = "primary",
                        solidHeader = TRUE,
                        plotOutput("plot_fd_overview", height = 300)
                    ),
                    box(
                        title = "Family Composition",
                        width = 6,
                        status = "info",
                        solidHeader = TRUE,
                        plotlyOutput("plot_family_composition", height = 300)
                    )
                ),
                
                fluidRow(
                    box(
                        title = "Community Weighted Means",
                        width = 12,
                        status = "success",
                        solidHeader = TRUE,
                        plotOutput("plot_cwm_overview", height = 300)
                    )
                )
            ),
            
            # TAB 2: Species Explorer
            tabItem(
                tabName = "species",
                h2("Species Explorer"),
                
                fluidRow(
                    box(
                        title = "Species Abundance",
                        width = 12,
                        status = "primary",
                        solidHeader = TRUE,
                        plotlyOutput("plot_species_abundance", height = 400)
                    )
                ),
                
                fluidRow(
                    box(
                        title = "Species Details",
                        width = 12,
                        status = "info",
                        solidHeader = TRUE,
                        DTOutput("table_species_details")
                    )
                )
            ),
            
            # TAB 3: Trait Analysis
            tabItem(
                tabName = "traits",
                h2("Trait Analysis"),
                
                fluidRow(
                    box(
                        title = "Trait Distribution",
                        width = 6,
                        status = "primary",
                        solidHeader = TRUE,
                        plotlyOutput("plot_trait_dist", height = 400)
                    ),
                    box(
                        title = "Trait vs Abundance",
                        width = 6,
                        status = "info",
                        solidHeader = TRUE,
                        plotlyOutput("plot_trait_abundance", height = 400)
                    )
                ),
                
                fluidRow(
                    box(
                        title = "Trait Correlations",
                        width = 12,
                        status = "success",
                        solidHeader = TRUE,
                        plotOutput("plot_trait_correlations", height = 500)
                    )
                )
            ),
            
            # TAB 4: Functional Diversity
            tabItem(
                tabName = "fd",
                h2("Functional Diversity"),
                
                fluidRow(
                    valueBoxOutput("box_fric"),
                    valueBoxOutput("box_feve"),
                    valueBoxOutput("box_fdis")
                ),
                
                fluidRow(
                    box(
                        title = "FD Metrics Comparison",
                        width = 12,
                        status = "primary",
                        solidHeader = TRUE,
                        plotOutput("plot_fd_metrics", height = 400)
                    )
                ),
                
                fluidRow(
                    box(
                        title = "Functional Space (PCA)",
                        width = 12,
                        status = "info",
                        solidHeader = TRUE,
                        plotlyOutput("plot_functional_space", height = 500)
                    )
                )
            ),
            
            # TAB 5: Community Comparison (multi-sample)
            tabItem(
                tabName = "compare",
                h2("Community Comparison"),
                
                if (length(sample_ids) > 1) {
                    tagList(
                        fluidRow(
                            box(
                                title = "FD Metrics Across Samples",
                                width = 12,
                                status = "primary",
                                solidHeader = TRUE,
                                plotOutput("plot_fd_comparison", height = 400)
                            )
                        ),
                        fluidRow(
                            box(
                                title = "CWM Comparison",
                                width = 12,
                                status = "info",
                                solidHeader = TRUE,
                                plotOutput("plot_cwm_comparison", height = 400)
                            )
                        )
                    )
                } else {
                    box(
                        width = 12,
                        status = "warning",
                        p("Community comparison requires multiple samples.",
                          br(),
                          "Only one sample detected in results.")
                    )
                }
            ),
            
            # TAB 6: Data Tables
            tabItem(
                tabName = "tables",
                h2("Data Tables"),
                
                fluidRow(
                    box(
                        title = "Complete Trait Data",
                        width = 12,
                        status = "primary",
                        solidHeader = TRUE,
                        DTOutput("table_complete_data")
                    )
                ),
                
                fluidRow(
                    box(
                        title = "Functional Diversity Metrics",
                        width = 6,
                        status = "info",
                        solidHeader = TRUE,
                        DTOutput("table_fd_metrics")
                    ),
                    box(
                        title = "Community Weighted Means",
                        width = 6,
                        status = "success",
                        solidHeader = TRUE,
                        DTOutput("table_cwm")
                    )
                )
            )
        )
    )
)

# =============================================================================
# SHINY APP SERVER
# =============================================================================

server <- function(input, output, session) {
    
    # Reactive data filtering
    sample_data <- reactive({
        if (length(sample_ids) > 1) {
            traits_combined %>% filter(sample_id == input$sample_select)
        } else {
            traits_combined
        }
    })
    
    sample_fd <- reactive({
        if (length(sample_ids) > 1) {
            fd_combined %>% filter(Sample == input$sample_select)
        } else {
            fd_combined
        }
    })
    
    # Value boxes - Overview
    output$box_n_species <- renderValueBox({
        valueBox(
            n_distinct(sample_data()$species_name),
            "Species Detected",
            icon = icon("dna"),
            color = "aqua"
        )
    })
    
    output$box_n_families <- renderValueBox({
        valueBox(
            n_distinct(sample_data()$Family),
            "Families",
            icon = icon("sitemap"),
            color = "green"
        )
    })
    
    output$box_total_abundance <- renderValueBox({
        valueBox(
            format(sum(sample_data()$abundance), big.mark = ","),
            "Total Reads",
            icon = icon("chart-line"),
            color = "yellow"
        )
    })
    
    # Value boxes - FD
    output$box_fric <- renderValueBox({
        valueBox(
            round(sample_fd()$FRic, 3),
            "FRic",
            icon = icon("cube"),
            color = "purple"
        )
    })
    
    output$box_feve <- renderValueBox({
        valueBox(
            round(sample_fd()$FEve, 3),
            "FEve",
            icon = icon("balance-scale"),
            color = "maroon"
        )
    })
    
    output$box_fdis <- renderValueBox({
        valueBox(
            round(sample_fd()$FDis, 3),
            "FDis",
            icon = icon("expand"),
            color = "olive"
        )
    })
    
    # Plots - Overview
    output$plot_fd_overview <- renderPlot({
        fd_long <- sample_fd() %>%
            select(FRic, FEve, FDiv, FDis, RaoQ) %>%
            pivot_longer(everything(), names_to = "Metric", values_to = "Value")
        
        ggplot(fd_long, aes(x = Metric, y = Value, fill = Metric)) +
            geom_col() +
            geom_text(aes(label = round(Value, 3)), vjust = -0.5) +
            theme_minimal() +
            scale_fill_brewer(palette = "Set2") +
            theme(legend.position = "none") +
            labs(y = "Value", x = "")
    })
    
    output$plot_family_composition <- renderPlotly({
        family_data <- sample_data() %>%
            group_by(Family) %>%
            summarise(Abundance = sum(abundance), .groups = "drop") %>%
            arrange(desc(Abundance)) %>%
            head(10)
        
        plot_ly(family_data, x = ~Family, y = ~Abundance, type = "bar",
                marker = list(color = "steelblue")) %>%
            layout(xaxis = list(title = ""), yaxis = list(title = "Abundance"))
    })
    
    output$plot_cwm_overview <- renderPlot({
        if (length(sample_ids) > 1) {
            sample_cwm <- cwm_combined %>% filter(Sample == input$sample_select)
        } else {
            sample_cwm <- cwm_combined
        }
        
        cwm_long <- sample_cwm %>%
            pivot_longer(-Sample, names_to = "Trait", values_to = "CWM")
        
        ggplot(cwm_long, aes(x = Trait, y = CWM, fill = Trait)) +
            geom_col() +
            theme_minimal() +
            scale_fill_viridis_d() +
            theme(
                legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            labs(x = "", y = "Community Weighted Mean")
    })
    
    # Species plots
    output$plot_species_abundance <- renderPlotly({
        species_abund <- sample_data() %>%
            group_by(species_name) %>%
            summarise(Abundance = sum(abundance), .groups = "drop") %>%
            arrange(desc(Abundance)) %>%
            head(20)
        
        plot_ly(species_abund, x = ~reorder(species_name, Abundance), 
                y = ~Abundance, type = "bar",
                marker = list(color = "coral")) %>%
            layout(
                xaxis = list(title = ""),
                yaxis = list(title = "Abundance"),
                margin = list(b = 100)
            )
    })
    
    # Trait plots
    output$plot_trait_dist <- renderPlotly({
        plot_ly(sample_data(), x = as.formula(paste0("~", input$trait_select)),
                type = "histogram", marker = list(color = "steelblue")) %>%
            layout(xaxis = list(title = input$trait_select),
                   yaxis = list(title = "Count"))
    })
    
    output$plot_trait_abundance <- renderPlotly({
        plot_data <- sample_data() %>%
            distinct(species_name, .keep_all = TRUE)
        
        plot_ly(plot_data, 
                x = as.formula(paste0("~", input$trait_select)),
                y = ~log10(abundance + 1),
                type = "scatter",
                mode = "markers",
                text = ~species_name,
                marker = list(size = 10, color = "darkgreen")) %>%
            layout(
                xaxis = list(title = input$trait_select),
                yaxis = list(title = "log10(Abundance + 1)")
            )
    })
    
    # Tables
    output$table_species_details <- renderDT({
        sample_data() %>%
            group_by(species_name, Family) %>%
            summarise(
                Abundance = sum(abundance),
                Spore_vol = first(Spore_vol_um3),
                Aspect = first(Aspect_ratio),
                Color = first(Color_score),
                .groups = "drop"
            ) %>%
            datatable(options = list(pageLength = 15))
    })
    
    output$table_complete_data <- renderDT({
        sample_data() %>%
            select(species_name, Family, abundance, 
                   Spore_vol_um3, Aspect_ratio, Color_score,
                   Wall_investment, Ornamentation_height_um) %>%
            datatable(options = list(pageLength = 20), filter = "top")
    })
    
    output$table_fd_metrics <- renderDT({
        datatable(fd_combined, options = list(pageLength = 10))
    })
    
    output$table_cwm <- renderDT({
        datatable(cwm_combined, options = list(pageLength = 10))
    })
}

# =============================================================================
# LAUNCH APP
# =============================================================================

cat("🚀 Launching Shiny dashboard...\n\n")
cat("Dashboard will open in your default web browser.\n")
cat("Press Ctrl+C to stop the server.\n\n")

shinyApp(ui = ui, server = server)
