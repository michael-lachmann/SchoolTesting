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
#            sliderInput("test.freq",
#                        "Test frequency:",
#                        min=2, max=90,value=30),
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




if( file.exists("res_list.Rda") & ( max(file.mtime("SchoolSim.R"),file.mtime("app.R"))<file.mtime("res_list.Rda") )  ) {load("res_list.Rda")} else {
res.list=list()
}
input=list(school.size=2000,test.freq=30, threshold=100, R0=1.5, inf.rate=20/2500*1e5/14, duration=90, school.size=200,symp.p=0.3)

# Define server logic required to draw a histogram
server <- function(input, output,clientData, session) {
    res.mass = reactive({
      
        n.days = 18
        n.types= 10
        N=as.numeric(input$school.size)
        x= setup.pop( n.days= n.days, n.types= n.types,symp.p = as.numeric(input$symp.p) )
        # How many infections will single individual produce?
        sum.inf = mean(colSums(x$inf[,x$inf.types])) # infection factor over days
        total.inf = N * sum.inf      # interaction with all kids in school
      
        test.freqs=c(2,3,7,14,30,90)
        
        duration=90
        
        ###### NEED TO INITIALIZE SCHOOL FOR EACH RUN
        # 
        # withProgress(message = 'Making plot', value = 0, {
        #     setProgress( 0)
#          clM<<-makeCluster(7)  
#          print("register")
#          registerDoParallel(clM)
        par.s=paste( input$school.size, input$test.cutoff, input$inf.rate)
        if( !(par.s %in% names(res.list)) ) {
        res=lapply( test.freqs, function(test.freq) {
          tt=Sys.time()
          N.runs=100
          res=sum( input$test.freq,  input$inf.rate, duration, as.numeric(input$test.cutoff)) # This is just so things recalculate
          res1=foreach( i=seq_len(N.runs), #.combine = cbind, .multicombine=T,
                       .export = c("run.pop","infect.pop","next.gen","test.pop","input","isolate","incProgress","N",
                                   "n.types","x","duration","total.inf"),
                       .packages = c("magrittr","roperators"),
                       .inorder=F
                      ) %do% {
            #isolate({
              x$pop["S",  ] =c( 0, rmultinom( 1, N,rep( 1, n.types)), 0)
              l=run.pop( N,x,max.T= duration, infect.rate = input$inf.rate/1e5, test.cutoff=as.numeric(input$test.cutoff),
                         mass.test.every = test.freq, test.every = 200,
                         par=list(beta = 0.6/ total.inf) # actual number is R0 divided by expected number of infections.
                         )
              c( N=N, 
                 Infect=N-tail(l$res,1)[,c("Sus")]-tail(l$res,1)[,c("Imp")], 
                 Imp=tail(l$res,1)[,c("Imp")], 
                 Q=tail(l$res,1)[,c("all.Q")],
                 MaxI=max(l$res[,"Inf"])
                 )
            #})
          }
          
          res2=foreach( i=seq_len(N.runs), #.combine = cbind, .multicombine=T,
                        .export = c("run.pop","infect.pop","next.gen","test.pop","input","isolate","incProgress","N",
                                    "n.types","x","duration","total.inf"),
                       .packages = c("magrittr","roperators"),
                       .inorder=F
          ) %do% {
            isolate({
              x$pop["S",  ] =c( 0, rmultinom( 1, N,rep( 1, n.types)), 0)
              l=run.pop( N,x,max.T= duration, infect.rate = input$inf.rate/1e5, test.cutoff=as.numeric(input$test.cutoff),
                         mass.test.every = test.freq, test.every = 200,
                         par=list(beta = 1.6/ total.inf) # actual number is R0 divided by expected number of infections.
              )
              c( N=N, 
                 Infect=N-tail(l$res,1)[,c("Sus")]-tail(l$res,1)[,c("Imp")], 
                 Imp=tail(l$res,1)[,c("Imp")], 
                 Q=tail(l$res,1)[,c("all.Q")],
                 MaxI=max(l$res[,"Inf"])
              )
            })
          }
          
          res3=foreach( i=seq_len(N.runs), #.combine = cbind, .multicombine=T,
                        .export = c("run.pop","infect.pop","next.gen","test.pop","input","isolate","incProgress","N",
                                    "n.types","x","duration","total.inf"),
                       .packages = c("magrittr","roperators"),
                       .inorder=F
          ) %do% {
            isolate({
              x$pop["S",  ] =c( 0, rmultinom( 1, N,rep( 1, n.types)), 0)
              l=run.pop( N,x,max.T= duration, infect.rate = input$inf.rate/1e5, test.cutoff=as.numeric(input$test.cutoff),
                         mass.test.every = test.freq, test.every = 200,
                         par=list(beta = 4/ total.inf) # actual number is R0 divided by expected number of infections.
              )
              c( N=N, 
                 Infect=N-tail(l$res,1)[,c("Sus")]-tail(l$res,1)[,c("Imp")], 
                 Imp=tail(l$res,1)[,c("Imp")], 
                 Q=tail(l$res,1)[,c("all.Q")],
                 MaxI=max(l$res[,"Inf"])
              )
            })
          }
          
 #           stopCluster(clM)
            res1.m=Reduce(rbind,res1) 
            res2.m=Reduce(rbind,res2) 
            res3.m=Reduce(rbind,res3) 
            res=abind( lo=res1.m,med=res2.m,hi=res3.m,along=3)
            
            x1   = res[,"Infect",] / res[,"Imp",]
            x2   = res[,"Q",]/(res[,"Infect",]+res[,"Imp",])
            x3 = (res[,"Infect",]+res[,"Imp",])/res[,"N",]
            
            res2=abind( res, IIrat=x1, Qrat=x2,InfRat=x3,along=2 )
            
            
            
            print(Sys.time()-tt)
            apply(res2,c(2,3),median)
        # } )
        })  # lapply test.freqs
        res= res %>% abind(along=3) %>% aperm(perm = c(3,1,2))
        dimnames(res)[[1]]=test.freqs
        res.list[[par.s]] <<- res
        save(res.list,file="res_list.Rda")
    }
      res.list[[par.s]]
    })
    output$distPlot <- renderPlot({

    })
    output$plot2 = renderPlot( height=1024,width=1024,res=150,{
      res=res.mass()

      test.freqs = dimnames(res)[[1]] %>% as.numeric()
      layout(cbind(1:2,3:4))
      plot(test.freqs,res[,"hi","Infect"]/res[,"hi","Imp"],type="l",ylab="In school infections per imported case",xlab="days between tests",log="x",xaxt="n")
      axis(1,test.freqs)
      lines(test.freqs,res[,"lo","Infect"]/res[,"lo","Imp"],col=2)
      lines(test.freqs,res[,"med","Infect"]/res[,"med","Imp"],col=3)
      legend("topleft",legend = c("Israel R0=4","Ireland R0=1.6","Germany R0=0.5"),lty=1,col=c(1,3,2))
      
      plot(test.freqs,res[,"hi","Q"]/(res[,"hi","Infect"]+res[,"hi","Imp"])*100,type="l",ylab="Percent cases caught",xlab="days between tests",xaxt="n",log="x",ylim=c(0,100))
      axis(1,test.freqs)
      lines(test.freqs,res[,"lo","Q"]/(res[,"lo","Infect"]+res[,"lo","Imp"])*100,col=2)
      lines(test.freqs,res[,"med","Q"]/(res[,"med","Infect"]+res[,"med","Imp"])*100,col=3)
      
      
      plot(test.freqs,res[,"hi","MaxI"],type="l",ylab="Maximal infected at once",xlab="days between tests",xaxt="n",log="x")
      axis(1,test.freqs)
      lines(test.freqs,res[,"lo","MaxI"],col=2)
      lines(test.freqs,res[,"med","MaxI"],col=3)
      
      plot(test.freqs,(res[,"hi","Infect"]+res[,"hi","Imp"])/res[,"hi","N"]*100,type="l",ylab="Percent infected",xlab="days between tests",xaxt="n",log="x",ylim=c(0,100))
      axis(1,test.freqs)
      lines(test.freqs,(res[,"lo","Infect"]+res[,"lo","Imp"])/res[,"lo","N"]*100,type="l",col=2)
      lines(test.freqs,(res[,"med","Infect"]+res[,"med","Imp"])/res[,"med","N"]*100,type="l",col=3)
    })
}


print("running shiny")
res=shinyApp(ui = ui, server = server)
print("after")
res