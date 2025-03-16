library(shiny)
library(bslib)
library(dplyr)
library(readr)
library(ranger)
library(leaflet)
library(sf)
library(geodata)

ui <- page_sidebar(
  title = "Species Distribution Model Generator",
  theme = bs_theme(bootswatch = "flatly", primary = "#3e8d63"),
  sidebar = sidebar(
    width = 300,
    fileInput("file", "Upload presence data (CSV)", accept = c(".csv")),
    selectInput("species_col", "Species column", choices = NULL),
    selectInput("lat_col", "Latitude column", choices = NULL),
    selectInput("lon_col", "Longitude column", choices = NULL),
    checkboxGroupInput(
      "bioclim_vars",
      "Bioclim Variables:",
      choices = c(
        "Annual Mean Temperature" = "bio_1",
        "Annual Precipitation" = "bio_12",
        "Max Temperature of Warmest Month" = "bio_5",
        "Min Temperature of Coldest Month" = "bio_6"
      ),
      selected = c("bio_1", "bio_12")
    ),
    numericInput("ntree", "Number of trees:", 500, min = 100, max = 2000),
    actionButton("run_model", "Run Model", class = "btn-success btn-lg w-100")
  ),
  layout_columns(
    card(card_header("Presence Data"), leafletOutput("map_presence", height = "350px")),
    card(card_header("Distribution Model"), leafletOutput("map_prediction", height = "350px"))
  ),
  layout_columns(
    card(card_header("Variable Importance"), plotOutput("variable_importance", height = "250px")),
    card(card_header("Model Info"), verbatimTextOutput("model_summary"))
  )
)

server <- function(input, output, session) {
  # Reactive values
  presence_data <- reactiveVal(NULL)
  sdm_model <- reactiveVal(NULL)
  prediction_raster <- reactiveVal(NULL)

  # Update selectors when file is uploaded
  observeEvent(input$file, {
    req(input$file)
    data <- read.csv(input$file$datapath)
    presence_data(data)

    updateSelectInput(session, "species_col", choices = names(data))
    updateSelectInput(session, "lat_col", choices = names(data))
    updateSelectInput(session, "lon_col", choices = names(data))

    # Auto-select columns based on common naming patterns
    for (col in names(data)) {
      if (grepl("species|taxa|name", tolower(col))) {
        updateSelectInput(session, "species_col", selected = col)
      }
      if (grepl("lat|latitude", tolower(col))) {
        updateSelectInput(session, "lat_col", selected = col)
      }
      if (grepl("lon|lng|longitude", tolower(col))) {
        updateSelectInput(session, "lon_col", selected = col)
      }
    }
  })

  # Display presence points map
  output$map_presence <- renderLeaflet({
    req(presence_data(), input$species_col, input$lat_col, input$lon_col)

    data <- presence_data()

    leaflet() %>%
      addProviderTiles(providers$Esri.WorldTopoMap) %>%
      addCircleMarkers(
        data = data,
        lng = ~lon,
        lat = ~lat,
        radius = 5,
        color = "#ff7f00",
        fillOpacity = 0.5,
        stroke = FALSE,
        popup = ~ paste(
          "<b>Species:</b>", species,
          "<br><b>Lat:</b>", lat,
          "<br><b>Lon:</b>", lon
        )
      )
  })

  # Run model when button is clicked
  observeEvent(input$run_model, {
    req(presence_data(), input$species_col, input$lat_col, input$lon_col, input$bioclim_vars)

    withProgress(message = "Building model...", value = 0, {
      # Process presence data
      incProgress(0.1, detail = "Processing presence data")
      data <- presence_data()

      # Create spatial points
      presence_pts <- data.frame(
        species = data[[input$species_col]],
        lon = data[[input$lon_col]],
        lat = data[[input$lat_col]]
      ) %>%
        st_as_sf(coords = c("lon", "lat"), crs = 4326)

      coords <- data %>%
        select(lon = !!input$lon_col, lat = !!input$lat_col) %>%
        na.omit()

      # Define study area with buffer
      min_lon <- min(coords$lon) - 1
      max_lon <- max(coords$lon) + 1
      min_lat <- min(coords$lat) - 1
      max_lat <- max(coords$lat) + 1


      study_extent <- ext(min_lon, max_lon, min_lat, max_lat)
      # Get environmental data
      incProgress(0.3, detail = "Downloading environmental data")
      tryCatch(
        {
          bioclim_list <- list()

          bio <- worldclim_global(var = "bio", res = 10, path = tempdir())

          bio <- bio[[c(1, 5, 6, 12)]]
          names(bio) <- stringr::str_remove_all(names(bio), "wc2.1_10m_")

          bio_filter <- bio[[c(input$bioclim_vars)]]
          env_stack <- crop(bio_filter, study_extent)

          # Generate pseudo-absence points
          incProgress(0.2, detail = "Generating background points")

          bg_points <- spatSample(env_stack,
            size = 1000,
            method = "random", na.rm = TRUE, as.points = TRUE
          )
          bg_coords <- crds(bg_points)
          colnames(bg_coords) <- c("lon", "lat")
          bg_env <- terra::extract(env_stack, bg_points)
          bg_env <- cbind(bg_coords, bg_env)
          bg_env$presence <- 0

          # Extract environmental data for presence points
          incProgress(0.1, detail = "Extracting environmental data")

          presence_env <- terra::extract(env_stack, coords)
          presence_env <- cbind(coords, presence_env)
          presence_env$presence <- 1

          # Combine for model training
          model_data <- rbind(presence_env, bg_env)
          model_data <- na.omit(model_data)

          # Build random forest model via ranger
          incProgress(0.2, detail = "Building random forest model")

          rf_model <- ranger::ranger(presence ~ .,
            data = model_data %>% select(presence, names(env_stack)),
            num.trees = 500, na.action = "na.omit", importance = "impurity"
          )

          sdm_model(rf_model)

          # Generate prediction
          incProgress(0.2, detail = "Generating prediction map")
          pred <- terra::predict(env_stack, rf_model, na.rm = T)

          prediction_raster(pred)

          # Show model summary
          output$model_summary <- renderPrint({
            rf_model
          })
        },
        error = function(e) {
          showNotification(paste("Error:", e$message), type = "error")
        }
      )
    })
  })

  # Display prediction map
  output$map_prediction <- renderLeaflet({
    req(prediction_raster())

    pred <- prediction_raster()
    data <- presence_data()

    # Color palette for predictions
    pal <- colorNumeric(
      palette = colorRampPalette(c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))(100),
      domain = c(0, 1)
    )

    leaflet() %>%
      addProviderTiles(providers$Esri.WorldTopoMap) %>%
      addRasterImage(pred, colors = pal, opacity = 0.7) %>%
      addCircleMarkers(
        data = data,
        lng = ~lon,
        lat = ~lat,
        radius = 3,
        color = "black",
        fillColor = "yellow",
        fillOpacity = 1,
        weight = 1
      ) %>%
      addLegend(
        position = "bottomright",
        pal = pal,
        values = c(0, 1),
        title = "Probability of Presence"
      )
  })

  # Variable importance plot
  output$variable_importance <- renderPlot({
    req(sdm_model())

    # imp <- importance(rf_model)

    imp <- ranger::importance(sdm_model())
    imp_df <- data.frame(
      Variable = names(imp),
      Importance = imp
    )

    # Create nicer labels for bioclim variables
    var_labels <- names(imp)

    imp_df$Label <- var_labels[imp_df$Variable]
    imp_df$Label[is.na(imp_df$Label)] <- imp_df$Variable[is.na(imp_df$Label)]

    # Sort by importance
    imp_df <- imp_df[order(imp_df$Importance, decreasing = TRUE), ]

    # Plot
    par(mar = c(8, 4, 2, 2))
    barplot(
      imp_df$Importance,
      names.arg = imp_df$Label,
      las = 2,
      cex.names = 0.8,
      col = "steelblue",
      main = "Variable Importance",
      ylab = "Mean Decrease in Gini"
    )
  })
}

shinyApp(ui, server)
