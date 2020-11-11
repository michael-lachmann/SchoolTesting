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
                        min=0, max=5,value=1.6,step=0.1)
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
server <- function(input, output) {
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
      
        
        withProgress(message = 'Making plot', value = 0, {
            setProgress( 0)
            sapply( 1:200,function(i) {
                if( i %% 20 == 0)
                    incProgress( 1/10, detail = paste("Running sims", i))

                l=run.pop( N,x,max.T= 100, infect.rate = input$inf.rate/1e5, 
                           mass.test.every = input$test.freq, close.thresh = input$threshold,test.every = 200,
                           par=list(beta = input$R0/ total.inf) # actual number is R0 divided by expected number of infections.
                           )
                c( dim( l$res)[1]*2, N-tail(l$res,1)[,3] )
            })
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




next.gen=function( pop, inf, det, weekend=F ) {
    gamma=1
    pop1=pop
    n.types = dim(pop)[2]-2
    n.days  = dim(pop)[1]-2
    # Move infected up by a day
    # deterministic
    i = ( 1:n.days)+1        
    sI = sum( pop * inf) # total infectivity - will effect each student in school
    mov = pop[ i,  ] 
    pop[ i  ,  ]= 0
    pop[ i+1,  ]= pop[ i+1, ] + mov

    if( sum( pop < 0) > 0) {
        print("problem")
#        browser()
    }
    # new infections within school only on weekdays
    if( ! weekend) {
      # new infections. s is total infection pressure
      i.S.ntypes=(1:n.types)+1
      new_e = rbinom(  n.types, pop[ "S", i.S.ntypes], (1-exp( log(1-sI))) )    # how many new infections we have in each type, among Sus. 
  
      pop[ "S" , i.S.ntypes] = pop["S",  i.S.ntypes] - new_e  # move out of S
      pop[ "I1", i.S.ntypes] = pop["I1", i.S.ntypes] + new_e  # add to infections on first day - I1
    }    
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


# Setup a population. 
# returns s list, where 
# pop - current state of population
#       row - time progression
#       state of individual - for example, can be different time courses of disease or different locations/rooms etc
# inf - infectiousness of each type
# symp  - 1 or 0 depending on whether sumptomatic or not
# 
# 
# input: symp.p - proportion of symptomatic individuals


setup.pop = function( n.days=10, n.types=50, symp.p=0.25 ) {
    # rows are progression of disease, columns are types of infection - how viral load increases
    pop= matrix(0,n.days+2,n.types+2)  # days: S, I1...IN, R  x types: N,T1...TM, Q
    dimnames(pop)=list( inf=c("S",paste0("I",1:(n.days)),"R"),type=c("N",paste0("T",1:(n.types)),"Q" ))
    
    # set up infectiveness of types - from Larremore
    inf.v = rep(0,dim(pop)[1])
    inf.v[1+4     ] = 0.6   # I4 is first infectious at 60%
    inf.v[1+(5:12)] = seq(1,0,len=8) # I5 to I11 declining infectiousness

    inf = matrix(inf.v, dim(pop)[1], dim(pop)[2])
    dimnames(inf)=dimnames(pop)
    inf[,"Q"]=0 # last type is quarantine
    inf[,"N"]=0 # first type hasn't decided yet
    
    # set symptoms
    symp = matrix(0, dim(pop)[1], dim(pop)[2])
    dimnames(symp)=dimnames(pop)
    symp[(1+5):(n.days+1),(1:(n.types*symp.p))+1]=1   # after I5 symptomatic, symp.p of types are symptomatic
    
    # detection level (viral load)
    det = matrix(0, dim(pop)[1], dim(pop)[2])
    det[(1:n.days)+1,(1:n.types)] = sapply( 1:n.types, function(i) {viral.load( n.days, rgamma(1,1.8)+0.2 , runif(1,7,11)   )}  )
    dimnames(det) = dimnames(pop)
    det[    ,"Q"] = 0  # viral load of quarantined is 0
    det[ "S",   ] = 0  # susceptible
    det[ "R",   ] = 0  # recovered
    
    list( pop=pop, det=det, inf=inf, symp=symp )
}


# infect.pop=function(pop, n.infect) {
#     for( i in seq_len(n.infect) ) {
#         if( pop["S","N"] > 0 ) {
#             infect.type = grep("^T",colnames(pop)) %>% sample(size=1)
#             pop["S","N"] = pop["S","N"]-1
#             pop["I1",infect.type] = pop["I1",infect.type]+1
#         }
#     }
#     pop
# }

infect.pop=function(pop, infect.rate) {
  n.types = dim(pop)[2]-2
  n.days  = dim(pop)[1]-2
  i.S.ntypes=(1:n.types)+1
  n.infect=rbinom( n.types, pop["S" ,i.S.ntypes], infect.rate)
  pop[ "S" , i.S.ntypes]  = pop["S" ,i.S.ntypes] -n.infect
  pop[ "I1", i.S.ntypes] = pop["I1",i.S.ntypes] +n.infect
  pop
}


run.pop = function(N,x,max.T=200, infect.rate=1/1000, mass.test.every=14, 
                   mass.test.fract=0.5, test.every=100, close.thresh=2, symp.thresh=0.5, use.weekend = T,
                   par=list(beta=0.)) {
    res     =matrix( 0, max.T, 5 ) # record results
    n.days  =dim(x$pop)[1]-2
    n.types =dim(x$pop)[2]-2
    
    
    for( t in 1:max.T ) {
        x$t = x$t+1
        
        # Infections from outside happen every day
        x$pop = infect.pop( x$pop, infect.rate)
        
        # Move infection forward and infect within school
        
        # Check if it is a weekend
        if( use.weekend & (t %% 7 < 5) )  # days 5 and 6 are weekend
          weekend = T
        else
          weekend = F
        
        x$pop = next.gen( x$pop, x$inf* par$beta, x$det, weekend = weekend )
        
        
        x = test.pop( x, x$symp, cutoff = symp.thresh, freq=1)
        #    N.quar = x$N.quar
        N.quar = sum( x$pop[(1:n.days)+1,"Q"])
        
        
        if( t %% mass.test.every==0) {
            x = test.pop( x, det = x$det, cutoff = 3, freq = 1/mass.test.fract)
            N.quar = N.quar + x$N.quar 
        }
        if( N.quar >= close.thresh) {
            res[t,] = c(sum( x$pop[ (1:n.days)+1, (1:n.types)+1 ] ), 
                        sum( x$pop[ (1:n.days)+1, "Q"]),
                        sum( x$pop[ "S"         , -m]),
                        sum( x$pop[             , "Q"]),
                        sum( x$pop[ "R"         , -m]))
            res=res[1:t,]
            break
        }
        
        m=n.types+2
        x = test.pop( x, x$det, cutoff = 5, freq = test.every)
        res[t,] = c(sum( x$pop[ (1:n.days)+1, (1:n.types)+1 ] ), 
                    sum( x$pop[ (1:n.days)+1, "Q"           ] ),
                    sum( x$pop[ "S"         , -m            ] ),
                    sum( x$pop[             , "Q"           ] ),
                    sum( x$pop[ "R"         , -m            ] )
        )
    }
    list( x=x, res=res)
}


# 
# 
# t.state = matrix( 0, n.types, n.time)
# 
# t.state.advance = function( x, s.in) {
#     n.dim(x)[2]
#     s.out =x[,n]
#     
# }
# 




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
