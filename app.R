#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(magrittr)
library(ggplot2)

library(doParallel)
library(parallel)

source("SchoolSim.R")

makeCluster.conditional=function( cl, no_cores = detectCores(logical = TRUE)-1) {
  #browser()
  cl.name = deparse(substitute(cl))
  if( !exists( cl.name,envir = parenv() ) ) {
    make.cl = T
  } else if( ! "cluster" %in% class(cl)) {
    make.cl=T
  } else make.cl=F
  
  if( make.cl) {
    makeCluster(no_cores)  
  } else {
    cl
  }
}

#cl = makeCluster.conditional( cl)
cl = makeCluster(3)
registerDoParallel(cl)

print("done")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("School testing"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("school.size",
                        "School size:",
                        min = 1,
                        max = 5000,
                        value = 2300),
            sliderInput("inf.rate",
                        "Infection import rate per 100k:",
                        min=1, max=400,value=40),
            sliderInput("threshold",
                        "Outbreak threshold:",
                        min=1, max=10, value=3),
            sliderInput("test.freq",
                        "Test frequency:",
                        min=2, max=30,value=14),
            sliderInput("R0",
                        "R0:",
                        min=0, max=5,value=1.6,step=0.1),
            radioButtons("country",label="R0 by country",
                         choices=c(
                           "South Korea 0.03"=0.03,
                           "Germany 0.6"=0.6,
                           "Israel 4.0"=4.0,
                           "Self chosen"=NA
                                   ))
        ),

        # Show a plot of the generated distribution
        mainPanel(
                plotOutput("distPlot"),
                br(),
                plotOutput("plot2")
            )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output,clientData, session) {
    print("here")
    onStop(function() {cat("Session stopped\n");stopCluster(cl)})
    observe( {
      updateTextInput(session, "R0", NULL, as.numeric(input$country))
    })
    res.mass = reactive({
        N=input$school.size
        n.days = 18
        n.types= 20
        x= setup.pop( n.days= n.days, n.types= n.types )
        # Each susceptible individual start in a certain type
        # Choose N individuals among the n.types
        x$pop["S",  ] =c( 0, rmultinom( 1, N,rep( 1, n.types)), 0)
        # How many infections will single individual produce?
        sum.inf = sum( x$inf[,"T1"]) # infection factor over days
        total.inf = N * sum.inf      # interaction with all kids in school
      
        
        ###### NEED TO INITIALIZE SCHOOL FOR EACH RUN
        
        withProgress(message = 'Making plot', value = 0, {
            setProgress( 0)
#          clM<<-makeCluster(7)  
#          print("register")
#          registerDoParallel(clM)

          tt=Sys.time()
          N.runs=500
          res=sum( input$test.freq, input$R0, input$inf.rate)
            res=foreach( i=seq_len(N.runs), #.combine = cbind, .multicombine=T,
                         .export = c("run.pop","infect.pop","next.gen","test.pop","input","isolate","incProgress"),
                         .packages = c("magrittr"),
                         .inorder=F
                        ) %dopar% {
              isolate({
#                if( i %% round(N.runs/20) == 0)
#                    incProgress( 1/20, detail = paste("Running sims", i))

                l=run.pop( N,x,max.T= 100, infect.rate = input$inf.rate/1e5, 
                           mass.test.every = input$test.freq, close.thresh = input$threshold,test.every = 200,
                           par=list(beta = input$R0/ total.inf) # actual number is R0 divided by expected number of infections.
                           )
                c( dim( l$res)[1]*2, N-tail(l$res,1)[,3] )
              })
            }
 #           stopCluster(clM)
            res=Reduce(cbind,res)
            print(Sys.time()-tt)
            res
        } )
    })
    output$distPlot <- renderPlot({
                
        hist( res.mass()[1,],br=seq(0,200,by=20),main="Days till outbreak detected",xlab="days")
        abline( v=median(res.mass()[1,]),col="red" )

    })
    output$plot2 = renderPlot({
        hist( res.mass()[ 2,],main="Number of students infected by time outbreak detected",xlab="students")
       abline( v=median( res.mass()[ 2,]),col="red" )
    })
}


print("running shiny")
res=shinyApp(ui = ui, server = server)
print("after")
res