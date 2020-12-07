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
library(gridExtra)
library(dplyr)
library(reshape)


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
if( F ) {
cl = makeCluster(3)
registerDoParallel(cl)
}

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
                        min=1, max=400,value=100),
            sliderInput("duration",
                        "Number of days (including weekend):",
                        min=1, max=150, value=90),
            sliderInput("test.freq",
                        "Test frequency:",
                        min=2, max=90,value=30),
            # sliderInput("R0",
            #             "R0:",
            #             min=0, max=5,value=1.6,step=0.1),
            radioButtons("test.cutoff",label="Test type",inline=T,
                          choices=c(
                            "PCR"=3,
                            "Pooled PCR"=4,
                            "Antigen"=5
                                    )),
            radioButtons("symp.p",label="Symptomatic freq",inline=T,selected="0.3",
                         choices=c(
                           "Lower (10%)"=0.1,
                           "Middle (30%)"=0.3,
                           "High & Adult (50%)"=0.5
                         ))
        ),

        # Show a plot of the generated distribution
        mainPanel(
#                plotOutput("distPlot"),
#                br(),
                plotOutput("plot2")
            )
    )
)


input=list(test.freq=30, threshold=100, R0=1.5, inf.rate=2/2500*1e5/14, duration=90, school.size=200)

# Define server logic required to draw a histogram
server <- function(input, output,clientData, session) {
    print("here")
if(F){    onStop(function() {cat("Session stopped\n");stopCluster(cl)})}
    observe( {
      updateTextInput(session, "R0", NULL, as.numeric(input$country))
    })
    res.mass = reactive({
        N=input$school.size
        n.days = 18
        n.types= 10
        x= setup.pop( n.days= n.days, n.types= n.types,symp.p = as.numeric(input$symp.p) )
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
          N.runs=50
          res=sum( input$test.freq,  input$inf.rate, input$duration, as.numeric(input$test.cutoff)) # This is just so things recalculate
          res1=foreach( i=seq_len(N.runs), #.combine = cbind, .multicombine=T,
                       .export = c("run.pop","infect.pop","next.gen","test.pop","input","isolate","incProgress"),
                       .packages = c("magrittr","roperators"),
                       .inorder=F
                      ) %do% {
            isolate({
              x$pop["S",  ] =c( 0, rmultinom( 1, N,rep( 1, n.types)), 0)
              l=run.pop( N,x,max.T= input$duration, infect.rate = input$inf.rate/1e5, test.cutoff=as.numeric(input$test.cutoff),
                         mass.test.every = input$test.freq, test.every = 200,
                         par=list(beta = 0.6/ total.inf) # actual number is R0 divided by expected number of infections.
                         )
              c( N=dim( l$res)[1]*2, Infect=N-tail(l$res,1)[,c("Sus")]-tail(l$res,1)[,c("Imp")], Imp=tail(l$res,1)[,c("Imp")] )
            })
          }
          
          res2=foreach( i=seq_len(N.runs), #.combine = cbind, .multicombine=T,
                       .export = c("run.pop","infect.pop","next.gen","test.pop","input","isolate","incProgress"),
                       .packages = c("magrittr","roperators"),
                       .inorder=F
          ) %do% {
            isolate({
              x$pop["S",  ] =c( 0, rmultinom( 1, N,rep( 1, n.types)), 0)
              l=run.pop( N,x,max.T= input$duration, infect.rate = input$inf.rate/1e5, test.cutoff=as.numeric(input$test.cutoff),
                         mass.test.every = input$test.freq, test.every = 200,
                         par=list(beta = 1.6/ total.inf) # actual number is R0 divided by expected number of infections.
              )
              c( N=dim( l$res)[1]*2, Infect=N-tail(l$res,1)[,c("Sus")]-tail(l$res,1)[,c("Imp")], Imp=tail(l$res,1)[,c("Imp")] )
            })
          }
          
          res3=foreach( i=seq_len(N.runs), #.combine = cbind, .multicombine=T,
                       .export = c("run.pop","infect.pop","next.gen","test.pop","input","isolate","incProgress"),
                       .packages = c("magrittr","roperators"),
                       .inorder=F
          ) %do% {
            isolate({
              x$pop["S",  ] =c( 0, rmultinom( 1, N,rep( 1, n.types)), 0)
              l=run.pop( N,x,max.T= input$duration, infect.rate = input$inf.rate/1e5, test.cutoff=as.numeric(input$test.cutoff),
                         mass.test.every = input$test.freq, test.every = 200,
                         par=list(beta = 4/ total.inf) # actual number is R0 divided by expected number of infections.
              )
              c( N=dim( l$res)[1]*2, Infect=N-tail(l$res,1)[,c("Sus")]-tail(l$res,1)[,c("Imp")], Imp=tail(l$res,1)[,c("Imp")] )
            })
          }
          
 #           stopCluster(clM)
            res1=Reduce(rbind,res1) %>% as.data.frame
            res2=Reduce(rbind,res2) %>% as.data.frame
            res3=Reduce(rbind,res3) %>% as.data.frame
            print(Sys.time()-tt)
            list( lo=res1,med=res2,hi=res3)
        } )
    })
    output$distPlot <- renderPlot({

    })
    output$plot2 = renderPlot( height=1024,width=1024,res=150,{
      res=res.mass()
      df.lo = melt(res$lo)
      df.med = melt(res$med)
      df.hi = melt(res$hi)
      x.lo= df.lo %>% filter(variable=="Infect" | variable=="Imp")
      x.med = df.med %>% filter(variable=="Infect" | variable=="Imp")
      x.hi = df.hi %>% filter(variable=="Infect" | variable=="Imp")
      
      x.max=max( x.lo$value, x.med$value, x.hi$value )
      
      bin = min( floor(x.max/10)+1, 5)
      
      g.lo = ggplot( df.lo %>% filter(variable=="Infect" | variable=="Imp"), aes(x=value,fill=variable)) + 
        labs( title= "Frequency of number infected, low spread")+
        ylab("Percentage of schools") +
        xlab("number infected") +
        coord_cartesian( xlim=c(0,x.max) ) +
        scale_fill_discrete(labels=c("Transmission in school","Imported infections")) +
        geom_histogram( aes(y=bin*..density..*100),binwidth=bin,position = "identity",alpha=0.5)+
        labs(fill = "Type of transmission")
      g.med = ggplot( df.med %>% filter(variable=="Infect" | variable=="Imp"), aes(x=value,fill=variable)) + 
        labs( title= "Frequency of number infected, medium spread")+
        xlab("number infected") +
        ylab("Percentage of schools") +
        coord_cartesian( xlim=c(0,x.max) ) +
        scale_fill_discrete(labels=c("Transmission in school","Imported infections")) +
        geom_histogram( aes(y=bin*..density..*100),binwidth=bin,position = "identity",alpha=0.5) +
        labs(fill = "Type of transmission")
      g.hi = ggplot( df.hi %>% filter(variable=="Infect" | variable=="Imp"), aes(x=value,fill=variable)) + 
        labs( title= "Frequency of number infected, high spread")+
        xlab("number infected") +
        ylab("Percentage of schools") +
        coord_cartesian( xlim=c(0,x.max) ) +
        scale_fill_discrete(labels=c("Transmission in school","Imported infections")) +
        geom_histogram( aes(y=bin*..density..*100),binwidth=bin,position = "identity",alpha=0.5) +
        labs(fill = "Type of transmission")
      grid.arrange(g.lo,g.med,g.hi,nrow=3)
    })
}


print("running shiny")
res=shinyApp(ui = ui, server = server)
print("after")
res