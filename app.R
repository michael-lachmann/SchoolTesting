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

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("School testing"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("School.size",
                        "School.size:",
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
                        min=2, max=30,value=14)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           plotOutput("plot2")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        N=input$School.size
        x=setup.pop()
        N = 2300
        x$pop["S","N"]=N
        withProgress(message = 'Making plot', value = 0, {
            
            res.mass=sapply(1:1000,function(i) {
                if( i %% 100==0)
                incProgress(1/10, detail = paste("Running sims", i))
                l=run.pop(N,x,max.T=100,infect.rate = input$inf.rate/1e5,mass.test.every = input$test.freq/2,close.thresh = input$threshold,test.every = 200)
                c( dim(l$res)[1]*2, N-tail(l$res,1)[,3] )
            })
        } )
                
        # res.rapid=sapply(1:1000,function(i) {
        #     l=run.pop(N,x,max.T=100,infect.rate = 40/1e5,mass.test.every = 200,close.thresh = 3,test.every = 1)
        #     c( dim(l$res)[1]*2, N-tail(l$res,1)[,3] )
        # })
        # 
        # res.symp=sapply(1:1000,function(i) {
        #     l=run.pop(N,x,max.T=100,infect.rate = 40/1e5,mass.test.every = 200,close.thresh = 3,test.every = 200)
        #     c( dim(l$res)[1]*2, N-tail(l$res,1)[,3] )
        # })
        # 
#        layout(rbind(1,2,3))
        hist(res.mass[1,],br=seq(0,200,by=20),main="Days till outbreak detected",xlab="days")
        abline( v=median(res.mass[1,]),col="red" )
        #hist(res.rapid[1,],br=seq(0,200,by=20))
        #hist(res.symp[1,],br=seq(0,200,by=20))
        
        
    })
    output$plot2 = renderPlot({
        hist(res.mass[2,],main="Number of students infected by time outbreak detected",xlab="students")
        abline( v=median(res.mass[2,]),col="red" )
    })
}




next.gen=function(pop, inf, det, gamma=0.5) {
    #    browser()  
    pop1=pop
    n.types = dim(pop)[2]-2
    n.days  = dim(pop)[1]-2
    # Move infected up by a day
    i = (1:n.days)+1        
    s = sum( pop * inf) /sum(pop)
    mov = rbinom(n.days*(n.types+2),pop[i, ],gamma) %>% array(.,dim=dim(pop[i,]))
    pop[i,  ] = pop[i,  ] - mov
    pop[i+1,] = pop[i+1,] + mov
    if( sum(pop<0) > 0) {
        print("problem")
#        browser()
    }
    # new infections. s is total infection pressure
    new_e = rbinom(  n.types, pop["S","N"], (1-exp(-s)) )    # how many new infections we have
    #    browser()
    new_e = rmultinom(1,new_e, rep(1/(n.types),n.types)) # distribute them among n.types
    pop["S","N"] = pop["S","N"] - sum(new_e)  # move out of S
    pop["I1",1+(1:n.types)] = pop["I1",1+(1:n.types)] + new_e  # add to infections on first day - I1
    
    # detect
    
    pop
}

test.pop=function( x, det, cutoff, freq, false.pos=0.00) {
    n.days  = dim( x$pop)[1]
    n.types = dim( x$pop)[2]
    i = cutoff < det
    quar = x$pop * 0 # empty copy of x$pop
    quar[i] = rbinom( length( x$pop[i]), x$pop[i], 1/freq) 
    x$pop[i] = x$pop[i] - quar[i] ;
    x$pop[, "Q"] = x$pop[, "Q"] + rowSums(quar)  # put in quarantine
    x$N.quar = sum(quar)
    quar = rbinom( 1, x$pop["S","N"], false.pos)
    x$pop["S","Q"] = x$pop["S","Q"]+quar
    x
}

viral.load=function(n,day.max, load.max) {
    day.max = min( day.max, n-1)
    x = c(seq( 0, load.max, by= load.max/day.max/10),
          seq( load.max, 0, by= load.max/(day.max-n)/10) )
    x[seq(1,length(x),len=n)]
}

ut=list()

input=list(beta1=0.1,beta2=0.1,gamma=1/7,tstart=1,socialDistStop=F)


N=1e6
n.types = 50
n.days = 10


setup.pop = function( n.days=10, n.type=50, symp.p=0.25 ) {
    # rows are progression of disease, columns are types of infection - how viral load increases
    pop= matrix(0,n.days+2,n.types+2)  # days: S, I1...IN, R  x types: N,T1...TM, Q
    dimnames(pop)=list( inf=c("S",paste0("I",1:(n.days)),"R"),type=c("N",paste0("T",1:(n.types)),"Q" ))
    
    # set us infectiveness of types
    inf = c(0,0.5,seq(1,0,len= n.days)) # 0 is for S, 0.5 for I1, and then goes down in the remaining n.days-1
    inf = matrix(inf, dim(pop)[1], dim(pop)[2])
    dimnames(inf)=dimnames(pop)
    inf[,"Q"]=0 # last type is quarantine
    inf[,"N"]=0 # first type hasn't decided yet
    
    # set symptoms
    symp = matrix(0, dim(pop)[1], dim(pop)[2])
    dimnames(symp)=dimnames(pop)
    symp[5:(n.days+1),(1:(n.type*symp.p))+1]=1   # after I4 symptomatic, 25% of cases are symptomatic
    
    # detection level (viral load)
    det = matrix(0, dim(pop)[1], dim(pop)[2])
    det[(1:n.days)+1,(1:n.types)] = sapply( 1:n.types, function(i) viral.load( n.days, rgamma(1,1.8)+0.2 , runif(1,7,11)   )  )
    dimnames(det) = dimnames(pop)
    det[,"Q"] = 0
    det["S",] = 0
    det["R",] = 0
    
    list( pop=pop, det=det, inf=inf, symp=symp )
}


infect.pop=function(pop, n.infect) {
    for( i in seq_len(n.infect) ) {
        if( pop["S","N"] > 0 ) {
            infect.type = grep("^T",colnames(pop)) %>% sample(size=1)
            pop["S","N"] = pop["S","N"]-1
            pop["I1",init.type] = pop["I1",infect.type]+1
        }
    }
    pop
}


run.pop = function(N,x,max.T=200, infect.rate=1/1000, mass.test.every=14, mass.test.fract=0.5, test.every=100, close.thresh=2, symp.thresh=0.5,
                   par=list(beta=0.28, gamma=1/2)) {
    #  browser()
    res=matrix(0,max.T,5)
    n.days=dim(x$pop)[1]-2
    n.types=dim(x$pop)[2]-2
    
    
    for( t in 1:max.T ) {
        x$t = x$t+1
        n.infect = rbinom(1,N, infect.rate/(7*par$gamma))
        x$pop = infect.pop( x$pop, n.infect)
        
        x$pop= next.gen( x$pop,x$inf* par$beta, x$det, par$gamma)
        
        x = test.pop( x, x$symp, cutoff = symp.thresh, freq=1)
        #    N.quar = x$N.quar
        N.quar = sum( x$pop[(1:n.days)+1,"Q"])
        
        if( t %% mass.test.every==0) {
            x = test.pop( x, det = x$det, cutoff = 3, freq = 1/mass.test.fract)
            N.quar = N.quar + x$N.quar 
        }
        if( N.quar >= close.thresh) {
            res[t,] = c(sum( x$pop[ (1:n.days)+1, (1:n.types)+1 ] ), 
                        sum(x$pop[(1:n.days)+1,"Q"]),
                        sum(x$pop["S",-m]),
                        sum(x$pop[,"Q"]),
                        sum(x$pop["R",-m]))
            res=res[1:t,]
            break
        }
        
        m=n.types+2
        x = test.pop( x, x$det, cutoff = 5, freq = test.every)
        res[t,] = c(sum( x$pop[ (1:n.days)+1, (1:n.types)+1 ] ), 
                    sum(x$pop[(1:n.days)+1,"Q"]),
                    sum(x$pop["S",-m]),
                    sum(x$pop[,"Q"]),
                    sum(x$pop["R",-m])
        )
    }
    list( x=x, res=res)
}

# x=setup.pop()
# N = 2300
# x$pop["S","N"]=N
# res.mass=sapply(1:1000,function(i) {
#     l=run.pop(N,x,max.T=100,infect.rate = 40/1e5,mass.test.every = 7,close.thresh = 3,test.every = 200)
#     c( dim(l$res)[1]*2, N-tail(l$res,1)[,3] )
# })
# 
# res.rapid=sapply(1:1000,function(i) {
#     l=run.pop(N,x,max.T=100,infect.rate = 40/1e5,mass.test.every = 200,close.thresh = 3,test.every = 1)
#     c( dim(l$res)[1]*2, N-tail(l$res,1)[,3] )
# })
# 
# res.symp=sapply(1:1000,function(i) {
#     l=run.pop(N,x,max.T=100,infect.rate = 40/1e5,mass.test.every = 200,close.thresh = 3,test.every = 200)
#     c( dim(l$res)[1]*2, N-tail(l$res,1)[,3] )
# })
# 
# layout(rbind(1,2,3))
# hist(res.mass[1,],br=seq(0,200,by=20))
# hist(res.rapid[1,],br=seq(0,200,by=20))
# hist(res.symp[1,],br=seq(0,200,by=20))
# 
# 
# layout(rbind(1,2,3))
# hist(res.mass[2,])
# hist(res.rapid[2,])
# hist(res.symp[2,])
# 

# Run the application 
shinyApp(ui = ui, server = server)
