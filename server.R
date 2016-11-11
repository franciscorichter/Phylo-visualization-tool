source('dependences.R')

##########################################

shinyServer(function(input, output) {

  randomVals <- eventReactive(input$goButton, {
    seed=round(runif(1,1,10000000))
    phyl2(tt=input$tt, lambda0=input$lambda,mu0=input$mu,K=input$K, seed=seed)

  })
  output$distPlot <- renderPlot({
    if (input$drop){
      dropex <- drop.fossil(randomVals()$newick) # drop extinct species
      plot(dropex, show.tip.label = input$tip, edge.color = input$color, edge.width = input$width, type=input$type, no.margin = input$margin)
    }
    else{
      plot(randomVals()$newick,show.tip.label = input$tip, edge.color = input$color, edge.width = input$width,type=input$type,no.margin=input$margin)
    }
    if (input$axis) axisPhylo()
  })
  output$ltt <- renderPlot({
    ltt(randomVals()$newick)
  })

  output$Text  <- renderText({
    p <- subplex(par = c(8,0.175,0.9),fn = llik,n = randomVals()$n, E = randomVals()$E, t = randomVals()$t)
    paste("The ML estimations for this tree would be",paste('lambda=',as.character(p$par[1])),paste('K=',as.character((p$par[1]-p$par[3])/p$par[2])),paste('mu=',as.character(p$par[3])),sep='\n')

  })

})
