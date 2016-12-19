
shinyUI(
  fluidPage(
    plotOutput("distPlot"),
    sidebarPanel(
      titlePanel('Simulation Parameters'),
      numericInput("lambda", "Lambda:", min = 0, max = 10, value = 0.8),
      numericInput("K", "K:", min = 1, max = 200, value = 40),
      numericInput("mu", "Mu:", min = 0, max = 2, value = 0.1),
      numericInput("tt", "Crown time:", min = 1, max = 30, value = 15),
      selectInput('model',"Model",c('Diversity-dependence', "Protracted")),
      br(),
      actionButton("goButton", "Simulate tree"),
      ("Click the button to simulate a phylogenetic tree under given parameters")
    ),
    sidebarPanel(
      titlePanel('Plot options'),
      checkboxInput("drop", "Drop extinct:", FALSE),
      checkboxInput("tip", "Show tip label:", FALSE),
      selectInput('color','Edge color',c('black','blue','red','green','purple','rainbow(7)')),
      numericInput('width','Edge width',value=1),
      selectInput('margin','Margin',c(F,T)),
      selectInput('axis','Include time axis',c(F,T)),
      selectInput('type','Type of phylogeny',c('phylogram','cladogram','fan','unrooted','radial'))

    ),
    sidebarPanel(
      verbatimTextOutput("Text")
    ),
    plotOutput('ltt')
  )
)
