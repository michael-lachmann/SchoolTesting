library(pryr)
library(roperators)

# Next
next.gen=function( x, weekend=F ) {
  gamma=1
  pop1=x$pop
  #browser()
  n.types = length( x$inf.types )
  n.days  = length( x$inf.days )
  # Move infected up by a day
  # deterministic
  i = ( 1:n.days)+1        
  sI = sum( x$pop * x$inf) # total infectivity - will effect each student in school
  mov = x$pop[ x$inf.days,  ] 
  x$pop[ x$inf.days  ,  ]= 0
  x$pop[ x$inf.days+1,  ]= x$pop[ x$inf.days+1, ] + mov
  
  if( sum( x$pop < 0) > 0) {
    print("problem")
    #        browser()
  }
  # new infections within school only on weekdays
  if( ! weekend) {
    # new infections. s is total infection pressure
    new_e = rbinom(  n.types, x$pop[ "S", x$inf.types], (1-exp( log(1-sI))) )    # how many new infections we have in each type, among Sus. 
    
#    x$pop[ "S" , x$inf.types] = x$pop["S",  x$inf.types] - new_e  # move out of S
#    x$pop[ "I1", x$inf.types] = x$pop["I1", x$inf.types] + new_e  # add to infections on first day - I1
     x$pop[ "S" , x$inf.types] %-=% new_e  # move out of S
     x$pop[ "I1", x$inf.types] %+=% new_e  # add to infections on first day - I1
  }    
  # detect
  
  x$pop
}


test.pop=function( x, det, cutoff, freq, false.pos=0.00) {
  n.days  = dim( x$pop)[1]
  n.types = dim( x$pop)[2]
  i = cutoff < det
  quar = x$pop * 0 # empty copy of x$pop
  quar[i] = rbinom( length( x$pop[i]), x$pop[i], 1/freq) 
  x$pop[i] %-=%  quar[i] ;
  x$pop[, "Q"] %+=%  rowSums(quar)  # put in quarantine
  x$N.quar = sum(quar)
  quar = rbinom( 1, x$pop["S","N"], false.pos)
  x$pop["S","Q"] %+=% quar
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

  inf.types= colnames(pop) %>% startsWith(prefix="T") %>% which
  inf.days = rownames(pop) %>% startsWith(prefix="I") %>% which

  
  # set up infectiveness of types - from Larremore
  inf.decline.t=n.days-5
  inf.v = rep(0,dim(pop)[1])
  inf.v[1+4     ] = 0.6   # I4 is first infectious at 60%
  inf.v[5+(1:inf.decline.t)] = seq(1,0,len=inf.decline.t) # I5 to I11 declining infectiousness
  
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
  det[ inf.days, inf.types] = sapply( 1:n.types, function(i) {viral.load( n.days, rgamma(1,1.8)+0.2 , runif(1,7,11)   )}  )
  dimnames(det) = dimnames(pop)
  det[    ,"Q"] = 0  # viral load of quarantined is 0
  det[ "S",   ] = 0  # susceptible
  det[ "R",   ] = 0  # recovered
#  as.environment( 
    list( pop=pop, det=det, inf=inf, symp=symp, inf.types=inf.types, inf.days=inf.days ) 
#  )
}


# infect.pop=function(pop, n.infect) 
#     for( i in seq_len(n.infect) ) {
#         if( pop["S","N"] > 0 ) {
#             infect.type = grep("^T",colnames(pop)) %>% sample(size=1)
#             pop["S","N"] = pop["S","N"]-1
#             pop["I1",infect.type] = pop["I1",infect.type]+1
#         }
#     }
#     pop
# }

infect.pop=function(x, infect.rate) {
  n.types = length(x$inf.types)
  n.days  = length(x$inf.days)
  n.infect=rbinom( n.types, x$pop[ "S" , x$inf.types], infect.rate)
  x$pop[ "S" , x$inf.types] %-=%  n.infect
  x$pop[ "I1", x$inf.types] %+=%  n.infect
  x$pop
}



run.pop = function(N,x,max.T=200, infect.rate=1/1000, mass.test.every=14, 
                   mass.test.fract=1, test.every=100, close.thresh=2, symp.thresh=0.5, use.weekend = T,
                   par=list(beta=0.)) {
  res     =matrix( 0, max.T, 5 ) # record results
  n.days  = length(x$inf.days)
  n.types = length(x$inf.types)
  x$inf = x$inf * par$beta
  
  for( t in 1:max.T ) {
    if( length(x$t)==0) x$t=0

    x$t = x$t+1
    
    # Infections from outside happen every day
    x$pop[] = infect.pop( x, infect.rate)
    
    # Move infection forward and infect within school
    
    # Check if it is a weekend
    if( use.weekend & (x$t %% 7 < 5) )  # days 5 and 6 are weekend
      weekend = T
    else
      weekend = F
    
    x$pop = next.gen( x, weekend = weekend )
    
    
    x = test.pop( x, x$symp, cutoff = symp.thresh, freq=1)
    #    N.quar = x$N.quar
    N.quar = sum( x$pop[(1:n.days)+1,"Q"])
    
    # Mass testing with PCR
    if( x$t %% mass.test.every==0) {
      x = test.pop( x, det = x$det, cutoff = 3, freq = 1/mass.test.fract)
      N.quar = N.quar + x$N.quar 
    }
    if( N.quar >= close.thresh) {
      res[t,] = c(sum( x$pop[ x$inf.days, x$inf.types ] ), 
                  sum( x$pop[ x$inf.days, "Q"]),
                  sum( x$pop[ "S"         , -m]),
                  sum( x$pop[             , "Q"]),
                  sum( x$pop[ "R"         , -m]))
      res=res[1:t,]
      break
    }
    
    m=n.types+2
    # individual testing with antigen
  #  x = test.pop( x, x$det, cutoff = 5, freq = test.every)
    res[t,] = c(sum( x$pop[ x$inf.days, x$inf.types ] ), # still infected
                sum( x$pop[ x$inf.days, "Q"           ] ), # in quaran
                sum( x$pop[ "S"         , -m            ] ), # sus
                sum( x$pop[             , "Q"           ] ), # all quaran
                sum( x$pop[ "R"         , -m            ] )  # recovered
    )
  }
  colnames(res)=c("Inf","in.Q","Sus","all.Q","Rec")
  list( x=x, res=res)
}



check.R0 = function( N, R0, n.days=18, n.type=20) {
  x= setup.pop( n.days= n.days, n.types= n.types )
  # Each susceptible individual start in a certain type
  # Choose N individuals among the n.types
  x$pop["S",  ] =c( 0, rmultinom( 1, N,rep( 1, n.types)), 0)
  # How many infections will single individual produce?
  sum.inf = sum( x$inf[,"T1"]) # infection factor over days
  total.inf = N * sum.inf      # interaction with all kids in school
  beta = R0/ total.inf
  
  p = x$pop
  
  p["I1","T1"] = 1
  x$inf = x$inf * par$beta
  res = sapply(seq_len(1000), function(j) {
  for( i in 1:n.days) {
    x$pop = next.gen( x, weekend = F )
    x$pop["I1","Q"] = sum(x$pop["I1", x$inf.types])
    x$pop["I1", x$inf.types] = 0
  }
  sum(x$pop[,"Q"])
  } )
  mean(res)
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