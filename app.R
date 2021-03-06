#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape)


library(doParallel)
library(parallel)
library(abind)

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
if( T ) {
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
          sliderTextInput("school.size","School size:", 
                  choices = c(25,50,100,500,2000,5000), 
                  selected = 2000,
                  animate = FALSE, grid = T, 
                  hide_min_max = FALSE, from_fixed = FALSE,
                  to_fixed = FALSE, from_min = NULL, from_max = NULL, to_min = NULL,
                  to_max = NULL, force_edges = FALSE, width = NULL, pre = NULL, 
                  post = NULL, dragRange = TRUE),
            # sliderInput("school.size",
            #             "School size:",
            #             min = 1,
            #             max = 5000,
            #             value = 2300),
            sliderInput("inf.rate",
                        "Infection import rate per 100k:",
                        min=1, max=400,value=100),
#            sliderInput("duration",
#                        "Number of days (including weekend):",
#                        min=1, max=150, value=90),
            sliderInput("test.freq",
                        "Test every:",
                        min=2, max=90,value=30),
             sliderInput("R0",
                         "R0:",
                         min=0, max=5,value=1.6,step=0.1),
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




if( file.exists("res_list.Rda") & ( max(file.mtime("SchoolSim.R"),file.mtime("app.R"))<file.mtime("res_list.Rda") )  ) {load("res_list.Rda")} else {
res.list=list()
}
input=list(school.size=2000,test.freq=30, threshold=100, R0=1.5, inf.rate=20/2500*1e5/14, duration=90, school.size=200,symp.p=0.3,
           test.cutoff=3)

# Define server logic required to draw a histogram
server <- function(input, output,clientData, session) {
    res.mass = reactive({
      par.s=paste( input$school.size, input$test.cutoff, input$inf.rate, input$symp.p)
      if( !is.null(res.list[[par.s]])) {
        res = res.list[[par.s]]
      } else {
   #     browser()
  #    res=sum( input$test.freq,  input$inf.rate,  as.numeric(input$test.cutoff), as.numeric(input$symp.p)) # This is just so things recalculate
      res=run.schools( school.size = as.numeric(input$school.size),
                       test.cutoff = as.numeric(input$test.cutoff),
                       inf.rate = as.numeric(input$inf.rate)/1e5,
                       symp.p = as.numeric(input$symp.p) )
      
      res.list[[par.s]] <<- res
       save(res.list,file="res_list.Rda")
      }
    res
            
    })
    
    res.single = reactive({
        res=run.single.school( school.size = as.numeric(input$school.size),
                         test.cutoff = as.numeric(input$test.cutoff),
                         inf.rate = as.numeric(input$inf.rate)/1e5,
                         R0= as.numeric(input$R0),
                         test.freq = as.numeric( input$test.freq ),
                         symp.p = as.numeric(input$symp.p) )
        
      
      res
      
    })
    
    output$distPlot <- renderPlot({

    })
    output$plot2 = renderPlot( height=1024,width=1024,res=70,{
      res=res.mass()
      dimnames(res)[[3]]=c("hi","med","lo")
#      res=aperm(res,c(1,3,2))
      test.freqs = dimnames(res)[[1]] %>% as.numeric()
      
      layout(cbind( c(1:2,5,5),c(3:4,5,5)))
#      plot(test.freqs,res[,"hi","Infect"]/res[,"hi","Imp"],type="l",ylab="In school infections per imported case",xlab="days between tests",log="x",xaxt="n")
      plot(test.freqs,res[,"IIrat", "hi"],type="l",ylab="In school infections per imported case",xlab="days between tests",log="x",xaxt="n")
      axis(1,test.freqs)
      lines(test.freqs,res[,"IIrat", "med"],col=2)
      lines(test.freqs,res[,"IIrat", "lo"],col=3)
#      lines(test.freqs,res[,"med","Infect"]/res[,"med","Imp"],col=2)
#      lines(test.freqs,res[,"lo","Infect"]/res[,"lo","Imp"],col=3)
      legend("topleft",legend = c("Israel R0=4","Ireland R0=1.6","Germany R0=0.5"),lty=1,col=c(1,2,3))
      
#      plot(test.freqs,res[,"hi","Q"]/(res[,"hi","Infect"]+res[,"hi","Imp"])*100,type="l",ylab="Percent cases caught",xlab="days between tests",xaxt="n",log="x",ylim=c(0,100))
      plot(test.freqs,res[,"Qrat", "hi"]*100,type="l",ylab="Percent cases caught",xlab="days between tests",xaxt="n",log="x",ylim=c(0,100))
      axis(1,test.freqs)
      lines(test.freqs,res[,"Qrat", "med"]*100,col=2)
      lines(test.freqs,res[,"Qrat", "lo"]*100,col=3)
#      lines(test.freqs,res[,"med","Q"]/(res[,"med","Infect"]+res[,"med","Imp"])*100,col=2)
#      lines(test.freqs,res[,"lo","Q"]/(res[,"lo","Infect"]+res[,"lo","Imp"])*100,col=3)
      
      
      plot(test.freqs,res[,"MaxI","hi"],type="l",ylab="Maximal infected at once",xlab="days between tests",xaxt="n",log="x")
      axis(1,test.freqs)
      lines(test.freqs,res[,"MaxI","med"],col=2)
      lines(test.freqs,res[, "MaxI","lo"],col=3)
      
#      plot(test.freqs,(res[,"hi","Infect"]+res[,"hi","Imp"])/res[,"hi","N"]*100,type="l",ylab="Percent infected",xlab="days between tests",xaxt="n",log="x",ylim=c(0,100))
      plot(test.freqs, res[,"InfRat","hi"]*100,type="l",ylab="Percent infected",xlab="days between tests",xaxt="n",log="x",ylim=c(0,100))
      axis(1,test.freqs)
      lines(test.freqs, res[,"InfRat","med"]*100,type="l",col=2)
      lines(test.freqs, res[,"InfRat","lo"]*100,type="l",col=3)
#      lines(test.freqs, (res[,"med","Infect"]+res[,"med","Imp"])/res[,"med","N"]*100,type="l",col=2)
#      lines(test.freqs, (res[,"lo","Infect"]+res[,"lo","Imp"])/res[,"lo","N"]*100,type="l",col=3)

      l=res.single()
      
      N=as.numeric(input$school.size)
      T=90
      plot( 1:90, N - l$res[,"Sus"],type="l")
      rect( 1:T-1, 0,                 1:T-0.2, l$res[,"new.Imp"]                   ,col="red")
      rect( 1:T-1, l$res[,"new.Imp"] ,1:T-0.2, l$res[,"new.Imp"]+l$res[,"new.Inf"]           )
      
    })
}


print("running shiny")
res=shinyApp(ui = ui, server = server)
print("after")
res