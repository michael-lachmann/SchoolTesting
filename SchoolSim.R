library(roperators)

# next.gen does 
# 1. disease progression
# 2. transmission
next.gen=function( x, do.infect=T ) {
  gamma=1
  pop1=x$pop

  
  n.types = length( x$inf.types )
  n.days  = length( x$inf.days )
  
  # Step 1: sum total infectivity
  sI = sum( x$pop * x$inf) # total infectivity - will effect each student in school
  
  # Step 2: deteministic progression of disease
  # Move infected up by a day
  # deterministic
  mov = x$pop[ x$inf.days,  ] 
  x$pop[ x$inf.days  ,  ]   = 0
  x$pop[ x$inf.days+1,  ] %+=%  mov
  
  
  # step 3: infect within the population
  if( do.infect) {
    x=infect.pop( x, 1- exp( log(1-sI)) )
  }
  x$pop
}


# external import of cases
infect.pop=function( x, infect.rate) {
  n.types = length(x$inf.types)
  n.days  = length(x$inf.days)
  
  # pop["S",] is a row of pre-determined typed. If infected, they stay in their column
  new_e = rbinom( n.types, x$pop[ "S" , x$inf.types], infect.rate)
  x$new.infect = sum( new_e ) 
  x$pop[ "S" , x$inf.types] %-=%  new_e
  x$pop[ "I1", x$inf.types] %+=%  new_e
  x
}

# Test the population
# det gives the chance for each type x day to be detected
# Those that are detected are put in isolation (Q)
# Note that det is not taken from x, because it can also be used by symptoms or
# other means. 
# cutoff determines cutoff to be detected,
# freq is chance to be tested.

test.pop=function( x, det, cutoff, freq, false.pos=0.00) {
  n.days  = dim( x$pop)[1]
  n.types = dim( x$pop)[2]
  i = cutoff < det
  quar = x$pop * 0 # make empty copy of x$pop
  # fill it with those there were detected as having covid
  quar[i] = rbinom( length( x$pop[i]), x$pop[i], 1/freq) 
  # remove from pop
  x$pop[i] %-=%  quar[i] ;
  # add to quarantine, on same day of disease progression
  x$pop[, "Q"] %+=%  rowSums(quar)  # put in quarantine
  #record
  x$new.quar = sum(quar)
  # Now add some false positives
  quar = rbinom( 1, x$pop["S","N"], false.pos)
  x$pop["S","Q"] %+=% quar
  x
}


run.pop = function(N,x,max.T=200, infect.rate=1/1000, mass.test.every=14, 
                   mass.test.fract=1, test.every=100, close.thresh=1e6, symp.thresh=0.5, use.weekend = T, test.cutoff=3,
                   par=list(beta=0.)) {
  res     =matrix( 0, max.T, 6 ) # record results
  n.days  = length(x$inf.days)
  n.types = length(x$inf.types)
  x$inf = x$inf * par$beta
  x$N.import = 0
  
  for( t in 1:max.T ) {
    if( length(x$t)==0) x$t=0
    
    x$t = x$t+1
    
    
    # Move infection forward and infect within school
    
    # Check if it is a weekend
    if( use.weekend & (x$t %% 7 < 5) )  # days 5 and 6 are weekend
      weekend = F
    else
      weekend = T
    
    x = infect.pop( x, infect.rate )
    x$N.import %+=% x$new.infect
    
    x$pop[] = next.gen( x, do.infect = !weekend )
    
    
    
    # put symptomatic cases in isolation
    x = test.pop( x, x$symp, cutoff = symp.thresh, freq=1)
    #    N.quar = x$N.quar
    N.quar = sum( x$pop[ x$inf.days,"Q"])
    
    # Mass testing with PCR
    if( x$t %% mass.test.every==0 ) {
      x = test.pop( x, det = x$det, cutoff = test.cutoff, freq = 1/mass.test.fract)
      N.quar = N.quar + x$new.quar 
    }
    
    if( T)  {
      if( N.quar >= close.thresh) {
        res[t,] = c(sum( x$pop[ x$inf.days, x$inf.types ] ), 
                    sum( x$pop[ x$inf.days, "Q"]),
                    sum( x$pop[ "S"         , -m]),
                    sum( x$pop[             , "Q"]),
                    sum( x$pop[ "R"         , -m]))
        res=res[1:t,]
        break
      }
    }
    
    m=n.types+2
    # individual testing with antigen
    #  x = test.pop( x, x$det, cutoff = 5, freq = test.every)
    res[t,] = c(sum( x$pop[ x$inf.days, x$inf.types ] ), # still infected
                sum( x$pop[ x$inf.days, "Q"           ] ), # in quaran
                sum( x$pop[ "S"         , -m            ] ), # sus
                sum( x$pop[             , "Q"           ] ), # all quaran
                sum( x$pop[ "R"         , -m            ] ),  # recovered
                x$N.import
               )
  }
  colnames(res)=c("Inf","in.Q","Sus","all.Q","Rec","Imp")
  list( x=x, res=res)
}

viral.load=function(n,day.max, load.max) {
  day.max = min( day.max, n-1)
  x = c(seq( 0, load.max, by= load.max/day.max/10),
        seq( load.max, 0, by= load.max/(day.max-n)/10) )
  x[seq(1,length(x),len=n)]
}

ut=list()

input=list(beta1=0.1,beta2=0.1,gamma=1/7,tstart=1,socialDistStop=F)


#N=1e6
#n.types = 50
#n.days = 10


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
  inf.v[1+3     ] = 0.6   # I4 is first infectious at 60%
  inf.v[1+4     ] = 1
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









check.R0 = function( N, R0, n.days=18, n.type=20) {
  x= setup.pop( n.days= n.days, n.types= n.types )
  # Each susceptible individual start in a certain type
  # Choose N individuals among the n.types
  x$pop["S",  ] =c( 0, rmultinom( 1, N,rep( 1, n.types)), 0)
  # How many infections will single individual produce?
  sum.inf = mean(colSums(x$inf[,x$inf.types])) # infection factor over days
  total.inf = N * sum.inf      # interaction with all kids in school
  beta = R0/ total.inf
  
  pop0=x$pop
  
  x$inf = x$inf * beta
  res = sapply(seq_len(1000), function(j) {
    x$pop= pop0
    x$pop["I1","T1"]=1
  for( i in 1:n.days) {
    x$pop = next.gen( x, do.infect = T )
#    x$pop["I1","Q"] = sum(x$pop["I1", x$inf.types])
#    x$pop["I1", x$inf.types] = 0
  }
  N-sum(x$pop["S",])-1
  } )
  mean(res)
}

  


run.schools=function( 
    n.days = 18,
    n.types= 10,
    school.size=2000,
    inf.rate=100/1e5,
    test.cutoff=3,
    symp.p=0.3,
    test.freqs=c(2,3,7,14,30,90),
    R0s = c(0.6,1.6,4),
    duration=90,
    N.runs=30
  ) {
  N= school.size
  x= setup.pop( n.days= n.days, n.types= n.types,symp.p = symp.p )
  
  # How many infections will single individual produce?
  sum.inf = mean(colSums(x$inf[,x$inf.types])) # infection factor over days
  total.inf = N * sum.inf      # interaction with all kids in school

  ###### NEED TO INITIALIZE SCHOOL FOR EACH RUN
  # 
  # withProgress(message = 'Making plot', value = 0, {
  #     setProgress( 0)
  #          clM<<-makeCluster(7)  
  #          print("register")
  #          registerDoParallel(clM)
  #        if(T) {
  res=lapply( test.freqs, function(test.freq) {
    res=lapply( R0s, function(R0) {
      foreach( i=seq_len(N.runs), #.combine = cbind, .multicombine=T,
                    .export = c("run.pop","infect.pop","next.gen","test.pop","input","isolate","incProgress","N",
                                "n.types","x","duration","total.inf"),
                    .packages = c("magrittr","roperators"),
                    .inorder=F
      ) %do% {
          x$pop["S",  ] =c( 0, rmultinom( 1, N,rep( 1, n.types)), 0)
          l=run.pop( N,x,max.T= duration, infect.rate = inf.rate, test.cutoff=test.cutoff,
                     mass.test.every = test.freq, test.every = 200,
                     par=list(beta = R0/ total.inf) # actual number is R0 divided by expected number of infections.
          )
          c( N=N, 
             Infect=N-tail(l$res,1)[,c("Sus")]-tail(l$res,1)[,c("Imp")], 
             Imp=tail(l$res,1)[,c("Imp")], 
             Q=tail(l$res,1)[,c("all.Q")],
             MaxI=max(l$res[,"Inf"])
          )
      }
    }) %>%
      lapply( FUN= function(x) Reduce(rbind,x) ) %>%  
      abind(  along=3)
    
    x1   = res[,"Infect",] / res[,"Imp",]
    x2   = res[,"Q",]/(res[,"Infect",]+res[,"Imp",])
    x3 = (res[,"Infect",]+res[,"Imp",])/res[,"N",]
    
    res2=abind( res, IIrat=x1, Qrat=x2,InfRat=x3,along=2 )
    apply(res2,c(2,3),median,na.rm=T)
    
  })  # lapply test.freqs
  res= res %>% abind(along=3) %>% aperm(perm = c(3,1,2))
  dimnames(res)[[1]]=test.freqs
  res
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