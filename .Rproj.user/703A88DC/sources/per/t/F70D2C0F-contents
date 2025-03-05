library(shiny)
library(bslib)
library(readr)
library(dplyr)
library(ggplot2)
library(terra)
library(leaflet)
library(randomForest)
library(elevatr)
library(geodata)

ui <- page_sidebar(
  title = "Species Distribution Modeling",
  sidebar = sidebar(
    fileInput("speciesFile", "Upload CSV file", accept = c(".csv")),
    selectInput("latColumn", "Latitude Column", choices = NULL),
    selectInput("lonColumn", "Longitude Column", choices = NULL),
    selectInput("speciesColumn", "Species Column", choices = NULL),
    numericInput("resolution", "Resolution (arc-minutes)", value = 10, min = 2.5, max = 30),
    checkboxGroupInput("bioclimVars", "Bioclimatic Variables",
                       choices = c("Annual Mean Temperature" = "tavg",
                                   "Annual Precipitation" = "prec",
                                   "Max Temperature of Warmest Month" = "tmax",
                                   "Min Temperature of Coldest Month" = "tmin"),
                       selected = c("tavg", "prec")),
    checkboxInput("includeElevation", "Include Elevation", TRUE),
    numericInput("ntrees", "Number of Trees", value = 500, min = 100, max = 2000),
    actionButton("runModel", "Run Model", class = "btn-primary")
  ),
  
  layout_columns(
    card(card_header("Data Preview"), tableOutput("dataPreview")),
    navset_card_tab(
      title = "Results",
      nav_panel("Map", leafletOutput("map", height = "500px")),
      nav_panel("Variable Importance", plotOutput("varImpPlot", height = "500px"))
    )
  )
)

server <- function(input, output, session) {
  species_data <- reactiveVal(NULL)
  rf_model <- reactiveVal(NULL)
  prediction_raster <- reactiveVal(NULL)
  
  observeEvent(input$speciesFile, {
    req(input$speciesFile)
    data <- read_csv(input$speciesFile$datapath)
    species_data(data)
    updateSelectInput(session, "latColumn", choices = colnames(data))
    updateSelectInput(session, "lonColumn", choices = colnames(data))
    updateSelectInput(session, "speciesColumn", choices = colnames(data))
  })
  
  output$dataPreview <- renderTable({
    req(species_data())
    head(species_data(), 10)
  })
  
  observeEvent(input$runModel, {
    req(species_data(), input$latColumn, input$lonColumn, input$speciesColumn, input$bioclimVars)
    
    data <- species_data()
    coords <- data %>%
      select(lon = !!input$lonColumn, lat = !!input$latColumn) %>%
      na.omit()
    
    # Determine extent with buffer
    min_lon <- min(coords$lon) - 1
    max_lon <- max(coords$lon) + 1
    min_lat <- min(coords$lat) - 1
    max_lat <- max(coords$lat) + 1
    
    withProgress(message = 'Processing...', value = 0, {
      # Download bioclimatic variables
      incProgress(0.2, detail = "Downloading climate data")
      bioclim_list <- list()
      for(var in input$bioclimVars) {
        bio <- worldclim_global(var = var, res = input$resolution, path = tempdir())
        bio_crop <- crop(bio, ext(min_lon, max_lon, min_lat, max_lat))
        bioclim_list[[var]] <- bio_crop
      }

      # Get elevation data if requested
      if(input$includeElevation) {
        incProgress(0.4, detail = "Downloading elevation data")
        elev <- get_elev_raster(locations = coords, z = input$resolution, src = "aws")
        elev <- rast(elev)
        names(elev) <- "elevation"
        bioclim_list[["elevation"]] <- elev
      }
      
      # Create environmental stack and extract at occurrence points
      incProgress(0.6, detail = "Preparing model data")
      env_stack <- rast(bioclim_list)
      
      occurrence_data <- terra::extract(env_stack, coords)
      occurrence_data <- cbind(coords, occurrence_data)
      occurrence_data$presence <- 1
      
      # Generate background points
      bg_points <- spatSample(env_stack, size = nrow(occurrence_data) * 3, 
                              method = "random", na.rm = TRUE, as.points = TRUE)
      bg_coords <- crds(bg_points)
      colnames(bg_coords) <- c("lon", "lat")
      bg_data <- terra::extract(env_stack, bg_points)
      bg_data <- cbind(bg_coords, bg_data)
      bg_data$presence <- 0
      
      # Combine data and train model
      model_data <- rbind(occurrence_data, bg_data)
      model_data <- na.omit(model_data)
      model_data$presence <- factor(model_data$presence)
      
      print(names(env_stack))
      
      incProgress(0.8, detail = "Training random forest model")
      model_formula <- as.formula(paste("presence ~", paste(names(env_stack), collapse = " + ")))
      
      rf <- randomForest(model_formula, data = model_data, 
                         ntree = input$ntrees, importance = TRUE)
      rf_model(rf)
      
      # Predict to the entire region
      prediction <- predict(env_stack, rf, type = "prob")
      prediction_raster(prediction)
      
      incProgress(1.0, detail = "Done!")
    })
  })
  
  output$map <- renderLeaflet({
    req(prediction_raster())
    
    pred_rast <- prediction_raster()
    
    # Get species data for plotting occurrence points
    data <- species_data()
    occurrence_points <- data %>%
      select(lon = !!input$lonColumn, lat = !!input$latColumn) %>%
      na.omit()
    
    # Create leaflet map
    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      addRasterImage(pred_rast, colors = colorRampPalette(c("blue", "green", "yellow", "red"))(100),
                     opacity = 0.7, group = "Prediction") %>%
      addCircleMarkers(data = occurrence_points, ~lon, ~lat, 
                       radius = 4, color = "black", fillColor = "white",
                       fillOpacity = 0.7, weight = 1.5, group = "Occurrences") %>%
      addLayersControl(
        overlayGroups = c("Prediction", "Occurrences"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  output$varImpPlot <- renderPlot({
    req(rf_model())
    
    var_imp <- importance(rf_model())
    var_imp_df <- data.frame(
      Variable = rownames(var_imp),
      MeanDecreaseAccuracy = var_imp[, "%IncMSE"]
    )
    
    var_imp_df <- var_imp_df[order(var_imp_df$MeanDecreaseAccuracy, decreasing = TRUE), ]
    
    # Make variable names more readable
    var_imp_df$Variable <- recode(var_imp_df$Variable,
                                  "tavg" = "Annual Mean Temperature",
                                  "tmax" = "Max Temp Warmest Month",
                                  "tmin" = "Min Temp Coldest Month",
                                  "prec" = "Annual Precipitation")
                                  # "elevation" = "Elevation")
    
    ggplot(var_imp_df, aes(x = reorder(Variable, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(x = "", y = "Variable Importance") +
      theme_minimal(base_size = 12)
  })
}

shinyApp(ui, server)