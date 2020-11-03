#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(expm)
library(plotly)
library(ggplot2)
library(ggplotify)
library(graphics)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Gene Tree Probability Calculator"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("length1",
                        "Branch length 1 in coalescent units:",
                        min = 0.001,
                        max = 10.0,
                        value = 1.0),
            sliderInput("length2",
                        "Branch length 2 in coalescent units:",
                        min = 0.001,
                        max = 10.0,
                        value = 1.0),
            br(),
            h4("The species tree for which gene tree probabilities will be computed is shown below. Use the sliders above to adjust the branch lenghs."),
            br(),
            img(src='SpeciesTree.png',height="60%", width="60%",align = "center"),
            br(),
            br(),
          
        ),

        # Show a plot of the generated distribution
        mainPanel(plotOutput(outputId="GeneTreeProbPlot"),
                  br(),
                  tags$h5("** = gene tree topology that matches that of the species tree, whether viewed as rooted or unrooted",style="color:steelblue"),
                  tags$h5("* = gene tree topology that matches that of the species tree when both are viewed as unrooted",style="color:maroon") 
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$GeneTreeProbPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        trees = c("** (a,(b,(c,d)))","* ((a,b),(c,d))","* (b,(a,(c,d)))","((a,c),(b,d))","((a,d),(b,c))","(a,(c,(b,d)))","(a,(d,(b,c)))","(b,(c,(a,d)))","(b,(d,(a,c)))","(c,(d,(a,b)))","(c,(a,(b,d)))","(c,(b,(a,d)))","(d,(a,(b,c)))","(d,(b,(a,c)))","(d,(c,(a,b)))")
       
        C3 = input$length1
        B3 = input$length2
        p1 = (1-exp(-B3))*(1-exp(-C3))+(1/3)*(1-exp(-B3))*exp(-C3)+exp(-B3)*(1/3)*(1-(3/2)*exp(-C3)+(1/2)*exp(-3*C3))+(1/3)*exp(-B3)*(1/3)*((3/2)*exp(-C3)-(3/2)*exp(-3*C3))+(1/18)*exp(-B3)*exp(-3*C3)
        p2 = (1/3)*(1-exp(-B3))*exp(-C3)+(1/3)*exp(-B3)*(1/3)*((3/2)*exp(-C3)-(3/2)*exp(-3*C3))+(2/18)*exp(-B3)*exp(-3*C3)
        p3 = (1/3)*(1-exp(-B3))*exp(-C3)+(1/3)*exp(-B3)*(1/3)*((3/2)*exp(-C3)-(3/2)*exp(-3*C3))+(1/18)*exp(-B3)*exp(-3*C3)
        p4 = (1/3)*exp(-B3)*(1/3)*((3/2)*exp(-C3)-(3/2)*exp(-3*C3))+(2/18)*exp(-B3)*exp(-3*C3)
        p5 = (1/3)*exp(-B3)*(1/3)*((3/2)*exp(-C3)-(3/2)*exp(-3*C3))+(2/18)*exp(-B3)*exp(-3*C3)
        p6 = exp(-B3)*(1/3)*(1-(3/2)*exp(-C3)+(1/2)*exp(-3*C3))+(1/3)*exp(-B3)*(1/3)*((3/2)*exp(-C3)-(3/2)*exp(-3*C3))+(1/18)*exp(-B3)*exp(-3*C3)
        p7 = exp(-B3)*(1/3)*(1-(3/2)*exp(-C3)+(1/2)*exp(-3*C3))+(1/3)*exp(-B3)*(1/3)*((3/2)*exp(-C3)-(3/2)*exp(-3*C3))+(1/18)*exp(-B3)*exp(-3*C3)
        p8 = (1/18)*exp(-B3)*exp(-3*C3)
        p9 = (1/18)*exp(-B3)*exp(-3*C3)
        p10 = (1/18)*exp(-B3)*exp(-3*C3)
        p11 = (1/3)*exp(-B3)*(1/3)*((3/2)*exp(-C3)-(3/2)*exp(-3*C3))+(1/18)*exp(-B3)*exp(-3*C3)
        p12 = (1/18)*exp(-B3)*exp(-3*C3)
        p13 = (1/3)*exp(-B3)*(1/3)*((3/2)*exp(-C3)-(3/2)*exp(-3*C3))+(1/18)*exp(-B3)*exp(-3*C3)
        p14 = (1/18)*exp(-B3)*exp(-3*C3)
        p15 = (1/18)*exp(-B3)*exp(-3*C3)
        
        
        probs=c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15)
        
        probdist = data.frame(cbind(trees,probs))
        probdist$probs <-as.numeric(probdist$probs)
        probdist$trees <-factor(probdist$trees,levels=c("** (a,(b,(c,d)))","* ((a,b),(c,d))","* (b,(a,(c,d)))","((a,c),(b,d))","((a,d),(b,c))","(a,(c,(b,d)))","(a,(d,(b,c)))","(b,(c,(a,d)))","(b,(d,(a,c)))","(c,(d,(a,b)))","(c,(a,(b,d)))","(c,(b,(a,d)))","(d,(a,(b,c)))","(d,(b,(a,c)))","(d,(c,(a,b)))"))
    
        p<-ggplot(data=probdist, aes(x=trees, y=probs)) +
            geom_col(fill="steelblue")+
            theme(axis.text.x = element_text(angle = 90,color=c("steelblue","maroon","maroon","black","black","black","black","black","black","black","black","black","black","black","black")),axis.text=element_text(size=18),
                  axis.title=element_text(size=18,face="bold"), title=element_text(size=18, face='bold'))+expand_limits(y=c(0,1)) +
            xlab("Gene Tree") + ylab("Probability") +
            ggtitle("Gene Tree Topology Probabilities for Selected Branch Lengths")
        p
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
