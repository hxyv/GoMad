# R shiny web app script
# Install library
if (!require("shiny")) install.packages("shiny")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("NGLVieweR")) install.packages("NGLVieweR")

# Load library
library(shiny)
library(ggplot2)
library(NGLVieweR)

# Define reusable function for compiling and running the Go program
run_go_program <- function(file_path, stepEM, timeNewton) {
  # Compile the Go program
  compile_result <- system("go build", intern = TRUE)
  if (length(compile_result) > 0) {
    message("Compilation Output: ", compile_result)
  }
  
  # Run the Go program
  run_result <- system(paste("./go", file_path, stepEM, timeNewton), intern = TRUE)
  if (length(run_result) > 0) {
    message("Runtime Output: ", run_result)
  }
}

ui <- fluidPage(
    theme = bs_theme(version = 5, bootswatch = "lux"),
  titlePanel(h1("GoMad: Molecular Dynamics Simulation", align = "center", style = "color: #2E86C1;")),
  sidebarLayout(
    sidebarPanel(
      h3("Input Options", style = "color: #34495E;"),
      fileInput("proteinFile", "Upload PDB File:", accept = c(".pdb")),
      p("Accepted file types: .pdb.", style = "font-size: 12px; color: #7B7D7D;"),
      textInput("fileURL", "Or Enter the URL of the PDB File:"),
      actionButton("inputProteinVisu", "Input Protein", class = "btn btn-primary"),
      hr(),
      h4("Simulation Parameters", style = "color: #34495E;"),
      numericInput("emStep", "Energy Minimization Steps:", value = 10, min = 2, step = 1),
      numericInput("simTime", "Simulation Time (fs):", value = 10, min = 2, step = 1),
      actionButton("runGoCode", "Run MD Simulation", class = "btn btn-primary")
    ),
    mainPanel(
        hr(),
      h3("Input Protein Structure", style = "color: #34495E;"),
      NGLVieweROutput("inputProtein", height = "500px"),
      hr(),
      h3("RMSD Result", style = "color: #34495E;"),
      plotOutput("outputPlot", height = "300px"),
      hr(),
      h3("Final Protein Structure", style = "color: #34495E;"),
      NGLVieweROutput("finalStructure", height = "500px"),
      hr(),
      h3("Superimposed Structure", style = "color: #34495E;"),
      NGLVieweROutput("superimposedStructure", height = "500px"),
      hr(),
      p("The above plot shows the RMSD of the protein over simulation time, and the structure is visualized using NGLVieweR.",
        style = "font-size: 12px; color: #7B7D7D;")
    ),
    position = "left"  # Sidebar on the left
  )
)
 
# creat the server
server <- function(input, output) {
    
    # Input protein visualization
    observeEvent(input$inputProteinVisu, {
        file_path <- NULL
        if (!is.null(input$proteinFile)) {
            file_path <- input$proteinFile$datapath
        } else if (input$fileURL != "") {
            tempFile <- tempfile(fileext = ".pdb")
            download.file(input$fileURL, tempFile, mode = "wb")
            file_path <- tempFile
        }
        req(file_path)
        
        output$inputProtein <- renderNGLVieweR({
            NGLVieweR(file_path) %>%
                addRepresentation("cartoon", param = list(color = "blue")) %>%
                setQuality("high") %>%
                stageParameters(backgroundColor = "white", zoomSpeed = 1)
        })
        
    })

  # Running MD
  observeEvent(input$runGoCode, {
    withProgress(message = "MD Simulation", value = 0, {
      incProgress(0.1, detail = "Initializing...")
      Sys.sleep(1)  
      
      file_path <- NULL
      if (!is.null(input$proteinFile)) {
        file_path <- input$proteinFile$datapath
      } else if (input$fileURL != "") {
        tempFile <- tempfile(fileext = ".pdb")
        download.file(input$fileURL, tempFile, mode = "wb")
        file_path <- tempFile
      }
      req(file_path)
      
      # Render the input structural visualization in Shiny
      output$inputProtein <- renderNGLVieweR({
          NGLVieweR(file_path) %>%
              addRepresentation("cartoon", param = list(color = "blue")) %>%
              setQuality("high") %>%
              stageParameters(backgroundColor = "white", zoomSpeed = 1)
      })
      
      stepEM <- input$emStep
      if (is.null(stepEM) || stepEM <= 0) {
          showNotification("Please provide a valid energy minimization steps (positive number).", type = "error")
          return()
      }
      
      timeNewton <- input$simTime
      if (is.null(timeNewton) || timeNewton <= 0) {
        showNotification("Please provide a valid simulation time (positive number).", type = "error")
        return()
      }
     
      
      incProgress(0.4, detail = "Running Go simulation...")
      
      
      Sys.sleep(2)  
      run_go_program(file_path, stepEM, timeNewton)
      
      incProgress(0.7, detail = "Loading RMSD data...")
      rmsd_data <- tryCatch({
        read.csv("result/RMSD.csv", header = FALSE)
      }, error = function(e) {
        showNotification("RMSD data could not be loaded. Check the Go program output.", type = "error")
        return(NULL)
      })
      req(rmsd_data)
      colnames(rmsd_data) <- c("Time", "RMSD")
      
      # Create RMSD plot
      incProgress(0.9, detail = "Generating plot...")
      rmsd_plot <- ggplot(rmsd_data, aes(x = Time, y = RMSD)) +
        geom_line(color = "#327fa8", size = 1) +
        labs(
          title = "RMSD Over Time",
          x = "Time (fs)",
          y = "RMSD (Å)"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)
        )
      
      # Save the plot as a PNG file
      ggsave("visualization/RMSD.png", plot = rmsd_plot, width = 10, height = 6, dpi = 300)
      
      # Render the plot in Shiny
      output$outputPlot <- renderPlot({
        rmsd_plot
      })
      
      # Render the final protein structural visualization in Shiny
      output$finalStructure <- renderNGLVieweR({
          NGLVieweR("result/output.pdb") %>%
              addRepresentation("cartoon", param = list(color = "orange")) %>%
              setQuality("high") %>%
              stageParameters(backgroundColor = "white", zoomSpeed = 1)
      })
      
      # Render the superimposed structural visualization in Shiny
      output$superimposedStructure <- renderNGLVieweR({
        NGLVieweR("result/output.pdb") %>%
              addRepresentation("cartoon", param = list(color = "orange")) %>%
              addStructure(file_path) %>%
              addRepresentation("cartoon", param = list(color = "blue")) %>%
              setSuperpose(
                  reference = 1, 
                  sele_reference = ":A", 
                  sele_target = ":A", 
                  superpose = TRUE
              ) %>% 
          setQuality("high") %>%
              stageParameters(backgroundColor = "white", zoomSpeed = 1)
      })
      
      incProgress(1, detail = "Complete!")
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)
