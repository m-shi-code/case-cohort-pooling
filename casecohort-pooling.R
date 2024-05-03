###########################
#### Functions used    ####
###########################

single.Rboot <- function(sim,ind,cov.lst='E'){
  sim <- sim[ind,]
  sim$subcohort <- sim$subcohort==1
  sim.cch <- sim[sim$subcohort|sim$status==1,]
  n.rep <- sim.cch$age.end-sim.cch$age.start
  sim.cch2 <- sim.cch[rep(seq_len(nrow(sim.cch)), n.rep),]
  sim.cch2$age <-as.numeric(unlist(sapply(1:nrow(sim.cch),function(v){(sim.cch$age.start[v]+1):sim.cch$age.end[v]})))
  sim.cch2$status1 <- 0
  sim.cch2$status1[cumsum(n.rep)] <- sim.cch$status
  sim.cch2 <- sim.cch2[sim.cch2$subcohort|sim.cch2$status1==1,]
  covs <- paste(cov.lst,collapse='+')
  summary(glm2<- glm(as.formula(paste("status1 ~as.factor(age)", covs,sep='+')),data=sim.cch2, family=binomial))

  boot.est<-summary(glm2)$coef[2,] 
  boot.est
  
}

pool.Rboot <- function(fam.id,ind,ccoh.pool,cov.lst='E'){
  fam.id<- fam.id[ind]
  
  row.id <- NULL
  for (id1 in fam.id) row.id <- c(row.id,which(ccoh.pool$pools==id1))
  ccoh.pool <- ccoh.pool[row.id,]
  
  ### create person-year format
  n.rep <- ccoh.pool$age.end-ccoh.pool$age.start
  ccoh.pool2 <- ccoh.pool[rep(seq_len(nrow(ccoh.pool)), n.rep),]
  ccoh.pool2$age <-as.numeric(unlist(sapply(1:nrow(ccoh.pool),function(v){(ccoh.pool$age.start[v]+1):ccoh.pool$age.end[v]})))
  ccoh.pool2$status1 <- 0
  ccoh.pool2$status1[cumsum(n.rep)] <- ccoh.pool$status
  ccoh.pool2$age[ccoh.pool2$p.size==1]<- ccoh.pool2$age[ccoh.pool2$p.size==1]+0.5
  ccoh.pool2<- ccoh.pool2[ccoh.pool2$subcohort|ccoh.pool2$status1==1,]
  
  dim(ccoh.pool2)
  
  covs <- paste(cov.lst,collapse='+')
  suppressWarnings(summary(glm3.frac1<- glm(as.formula(paste("status1 ~as.factor(age)", covs,sep='+')),data=ccoh.pool2, weights=frac, family=binomial)))
  coef(glm3.frac1)[cov.lst]
}


#' @param study_week study week number
#' @param pool.size  pool size
#' @param n.cpu number of CPU used in bootstrapping
#' @param n.data number of simulated datasets
#' @param subcohort.frac fraction of the subcohort of the entire cohort
#' @param base.risk baseline risk of disease
#' @param beta.e log odds ratio of exposure


###############################
#### 1: Required Packages  ####
##############################

  # install required packages
  requiredPackages = c('survival','sandwich','lmtest','boot')
  for (pack1 in requiredPackages) {
    if (!require(pack1, character.only = TRUE))
      install.packages(pack1)
    library(pack1, character.only = TRUE)
  }
  

###################################
#### 2: Case-cohort simulation ####
###################################

  ## case-cohort simulation
  set.seed(177331)
  pool.size <- 4  ### pool size
  n.cpu <- 5  ### number of CPU used in bootstrapping
  n.data <- 40000 ### Number of simulated datasets
  subcohort.frac <- 0.04 ### fraction of the subcohort of the entire cohort.
  n.boot <- 10  ### number of bootstraps. Set to 10 for testing but recommend 1000
  base.risk <- 1/500  ###baseline risk
  beta.e <- log(1.05)  ###log odds ratio of exposure

  #### The following are parameters for the breast cancer simulation scenario
  beta.age <- log(1.01) ## log OR of age
  age.mid <- 50 ## risk decrease for age<50 and increase for age>50
  age.max <- 80 ## max age of followup
  n.fu <- 15  ## number of years if followup
  dropoff.rate <- 0.02 ## dropout rate
  
  
  ID <- 1:n.data
  age.start <- sample(35:74,n.data,replace=T)
  E<- exp(rnorm(n.data*1.2))
  E <- E[E<= log(1/base.risk)/beta.e]
  E <- E[1:n.data]
  
  sim<- data.frame(ID,age.start,age.end=age.start, E, status=rep(0,n.data))
  sim <- sim[sample(1:n.data),]
  n.fu <- 15
  
  id.drop <- NULL
  sim.all <- NULL
  for (i in 1:n.fu){
    probs <- exp(beta.e*sim$E+(sim$age.start+i-age.mid)*beta.age)*base.risk/(1+exp(beta.e*sim$E+(sim$age.start+i-age.mid)*beta.age)*base.risk)
    risk1<- lapply(probs, function(v){return(sample(c(0,1),1,prob=c(1-v,v)))})
    
    sim$age.end <- sim$age.end+1
    sim$status[which(risk1==1)]<-1
    sim.all <- rbind(sim.all, sim[sim$age.end==age.max | sim$status==1,] )
    sim <- sim[-which(sim$age.end==age.max | sim$status==1),]
    ## randomly pick  to drop out of study
    dropped <- sample(1:nrow(sim),dropoff.rate*nrow(sim))
    sim.all <- rbind(sim.all, sim[dropped,])
    id.drop<- c(id.drop, sim$ID[dropped])
    sim <- sim[-dropped, ]
  }
  sim<- rbind(sim, sim.all) 
  sim <- sim[order(sim$ID),]
  sim$n.fu <- sim$age.end-sim$age.start
  
  subcohort <- sample(sim$ID,n.data*subcohort.frac)
  ## sim is the entire cohort 
  sim$subcohort <- sim$ID%in%subcohort
  ## ccoh.data is the case-subcohort
  ccoh.data <- sim[sim$subcohort|sim$status==1,]

  ### ccoh.data variables: 
  ## ID: sample ID
  ## age.start: age at enrollment
  ## age.end: age at the end of followup for non-cases in the subcohort and age at diagnosis for all cases
  ## E: exposure value
  ## status: disease status at the end of follow up
  ## n.fu: years of follow up
  ## subcohort: indicator variable for subcohort status
  ## case: 1=non-subocohort cases
  
  ###################################################
  #### 2.1: Case-cohort analysis without pooling ####
  ###################################################
  
  ########################################
  #### 2.1.1: Create year by year data  ####
  ########################################
  
    n.rep <- sim$age.end-sim$age.start
    sim2 <- sim[rep(seq_len(nrow(sim)), n.rep),]
    sim2$age <-as.numeric(unlist(sapply(1:nrow(sim),function(v){(sim$age.start[v]+1):sim$age.end[v]})))
    sim2$status1 <- 0
    sim2$status1[cumsum(n.rep)] <- sim$status
    ccoh.data2 <- sim2[sim2$subcohort|sim2$status1==1,]  ## nonsubcohort cases only contribute the year at dx
    
  ##############################################
  #### 2.1.2: Logistic regression to analyze  ####
  ##############################################
    
   ###### entire cohort analysis 
   summary(glm1<- glm(status1 ~E+as.factor(age),data=sim2, family=binomial))
  
   ###### case-cohort analysis
   summary(glm2<- glm(status1 ~E+as.factor(age),data=ccoh.data2, family=binomial))
   
   ## bootstrap 
   beta.boot <- boot(sim,single.Rboot,R = n.boot, stype = "i", parallel = "multicore",ncpus=n.cpu)$t
   est.boot1 <- c(mean(beta.boot[,1]),sd(beta.boot[,1]),quantile(beta.boot[,1],c(0.025,0.975)))
   
   
   ############################################
   #### 3: Case-cohort pooling strategy    ####
   ############################################
   
   #######################################
   #### non-subcohort case pooling: pool based on age at diagnosis 
   #######################################
   
    ccoh.data$case <- 1*(ccoh.data$status==1 & !ccoh.data$subcohort) ## non-subcohort cases
    ccoh.case <- ccoh.data[ccoh.data$case==1,]
    ccoh.case <- ccoh.case[order(ccoh.case$age.end),]
    n.age.end <- table( ccoh.case$age.end)
    
    pool.ID <-0 ## keep track of last pool ID pools- IDs of all pools
    
    pools <- NULL #list of IDs either pooled or singles
    p.size <- NULL
    for (a in 1:length(n.age.end)){
      n.pool1 <- floor(n.age.end[a]/pool.size) # pool as much as possible
      id1<- NULL
      if (n.pool1>0) id1<- rep(1:n.pool1+pool.ID,each=pool.size)
      pool.ID <- n.pool1+pool.ID
      n.single1 <- n.age.end[a]-n.pool1*pool.size
      if (n.single1>0) id1<- c(id1,(1:n.single1)+pool.ID)
      pool.ID <- n.single1+pool.ID
      pools <- c(pools,id1)
      p.size <- c(p.size, c(rep(pool.size,n.pool1*pool.size),rep(1,n.single1)))
    }
    
    ccoh.case$pools <- pools
    ccoh.case$p.size <- p.size
    
    
    ######################################
    ##### Subcohort pooling: pool based on age at entry- cases within subcort will be followed up after dx ###
    ######################################
   ccoh.subc <-  ccoh.data[ccoh.data$case==0,]
   ccoh.subc <- ccoh.subc[order(ccoh.subc$age.start),]
   n.age.end <- table(ccoh.subc$age.start)  ## n.age.end here is actually n.age.start
    
    pools <- NULL
    p.size <- NULL
    for (a in 1:length(n.age.end)){
      n.pool1 <- floor(n.age.end[a]/pool.size) 
      id1<- NULL
      if (n.pool1>0) id1<- rep(1:n.pool1+pool.ID,each=pool.size)
      pool.ID <- n.pool1+pool.ID
      n.single1 <- n.age.end[a]-n.pool1*pool.size
      if (n.single1>0) id1<- c(id1,(1:n.single1)+pool.ID)
      pool.ID <- n.single1+pool.ID
      pools <- c(pools,id1)
      p.size <- c(p.size, c(rep(pool.size,n.pool1*pool.size),rep(1,n.single1)))
    }
    
    ccoh.subc$pools <- pools
    ccoh.subc$p.size <- p.size
    
    pool.mix <-tapply(ccoh.subc$status, ccoh.subc$pools, function(v) length(table(v))) ## check how many pools have cases in the subcohort pools
    ## pool IDs with both cases and noncases in the same pool
    pool.mixed.id <- as.numeric(names(pool.mix))[which(pool.mix>1)]

    ## cases in subcort age  at last fu minus one, except for the cases in a single pool
    ccoh.subc$age.end[ccoh.subc$status==1 & ccoh.subc$p.size!=1]<- ccoh.subc$age.end[ccoh.subc$status==1 & ccoh.subc$p.size!=1]-1
    status.subcohort <- rep(0,length(table(ccoh.subc$pools)))
    status.subcohort[unique(ccoh.subc$pools)%in%ccoh.subc$pools[ccoh.subc$status==1 &ccoh.subc$p.size==1]] <-1
    
    ccoh.pool <- as.data.frame(rbind(ccoh.case,ccoh.subc))
    write.table(ccoh.pool[,c('ID','pools','p.size')],row.names = F,sep='\t',quote=F,file='pool_mapping.txt')
    
 

    ###########################################
    #### 4: Form Case-cohort pools          ####
    ###########################################
    
    
    ## determine start and end age of each pool, create partial pools ##
    
    E<- tapply(ccoh.pool$E, ccoh.pool$pools, sum)
    E.pool <- data.frame(pools=names(E),E=E)
    write.table(E.pool,row.names = F,sep='\t',quote=F,file='Exposure_in_pools.txt')
    ### in real study pooled exposure values will be measured and saved in "Exposure_in_pools.txt" file
    ### read in exposure value and merge the exposure data with the rest of the variables
   
    age.start <- tapply(ccoh.pool$age.start,ccoh.pool$pools,max)
    age.end <- tapply(ccoh.pool$age.end,ccoh.pool$pools,min)
    
    status <- c(rep(1,length(table(ccoh.case$pools))),status.subcohort)
    subcohort <-  ccoh.pool$subcohort[!duplicated(ccoh.pool$pools)]
    pools <-  ccoh.pool$pools[!duplicated(ccoh.pool$pools)]
    p.size <-  ccoh.pool$p.size[!duplicated(ccoh.pool$pools)]
    
    ccoh.pool <- data.frame(E, age.start,age.end, status, subcohort,pools,p.size)
    
    ccoh.pool <- ccoh.pool[ccoh.pool$age.end> ccoh.pool$age.start,]
    ccoh.pool$frac <- 1
    remove(E, age.start,age.end, status, subcohort,pools,p.size)

    ####################################################################
    ## add fractional numbers for the subcohort pools with different end age due to dropout
    ####################################################################
    
    ## ccoh.subc is the subcohort individual samples
    
    ccoh.pool.org <- as.data.frame(ccoh.subc)
    n.end.age <- tapply(ccoh.pool.org$age.end,ccoh.pool.org$pools,function(v){length(table(v))})
    sum(n.end.age>1 )
    id.mixed.end <- names(n.end.age[n.end.age>1]) ## pools with different end age
    length(pool.mixed.id)
    length(id.mixed.end)

    partial.pool <- NULL
    
    for (id in id.mixed.end){
      pool1 <- ccoh.pool.org[ccoh.pool.org$pools==id,]
      E1 <- ccoh.pool$E[ccoh.pool$pools==id]
      age.start1 <- pool1$age.start[1]
      status1 <- pool1$status[1]
      subcohort1 <- pool1$subcohort[1]
      pools1 <- pool1$pools[1]
      p.size1 <- pool1$p.size[1]
      ages <- sort(unique(pool1$age.end))
      if (length(ages)>1){
        for (j in 2:length(ages)){
          age1 <- ages[j-1]
          age2 <- ages[j]
          
          frac1 <- sum(pool1$age.end >= age2)/nrow(pool1)
          if (frac1<1 & sum(pool1$status[pool1$age.end<age2])==0) partial.pool <- rbind(partial.pool,c(E1,age1,age2,status1,subcohort1,pools1,p.size1,frac1)) 
        }
      }
    }
    
    ccoh.pool$frac <-1
  if (!is.null(partial.pool)) {
    colnames(partial.pool) <- colnames(ccoh.pool)
    ccoh.pool <- rbind(ccoh.pool, partial.pool)
    }
    ccoh.pool <- ccoh.pool[order(ccoh.pool$pools),]
    
    
    ####################################################################
    #### 5. create person-year format                               ####
    ####################################################################
    
    n.rep <- ccoh.pool$age.end-ccoh.pool$age.start
    ccoh.pool2 <- ccoh.pool[rep(seq_len(nrow(ccoh.pool)), n.rep),]
    ccoh.pool2$age <-as.numeric(unlist(sapply(1:nrow(ccoh.pool),function(v){(ccoh.pool$age.start[v]+1):ccoh.pool$age.end[v]})))
    ccoh.pool2$status1 <- 0
    ccoh.pool2$status1[cumsum(n.rep)] <- ccoh.pool$status
    ### All singletons add half year to age so that they will have different stratification parameters
    ccoh.pool2$age[ccoh.pool2$p.size==1]<- ccoh.pool2$age[ccoh.pool2$p.size==1]+0.5
    ccoh.pool2<- ccoh.pool2[ccoh.pool2$subcohort|ccoh.pool2$status1==1,]
 

    
    #####################################################################
    #### 6. Logistic regression to analyze pooled case-cohort data   ####
    #####################################################################
    
    summary(glm3.frac1<- glm(status1 ~E+as.factor(age),data=ccoh.pool2, weights=frac, family=binomial))
    coef1 <- coeftest(glm3.frac1, vcov. = vcovCL(glm3.frac1, cluster = ccoh.pool2$pools, type = "HC0"))['E',]
    

    ####################################################################
    #### 7. Bootstrap method to get variance                        ####
    ####################################################################
    
    boot1 <- boot(unique(ccoh.pool$pools),pool.Rboot,R = n.boot, stype = "i", parallel = "multicore",ncpus=n.cpu,ccoh.pool=ccoh.pool)
    beta.boot.pool <- boot1$t
    est.boot1 <- c(mean(beta.boot.pool[,1]),sd(beta.boot.pool[,1]),quantile(beta.boot.pool[,1],c(0.025,0.975)))
    
    OR= paste(round(exp(coef1[1]),2)," (",round(exp(coef1[1]-1.96*est.boot1[2]),2),", ",round(exp(coef1[1]+1.96*est.boot1[2]),2),")",sep='')
    OR


    ###########################################################################################
    #
    #    Functions for real data analysis
    #
    ###########################################################################################

    ############################################################
    #### Function to form pools                             ####
    ############################################################
create.pools <- function(ccoh.data,cov.lst=NA){

        if (sum(!(c('ID', 'age.start', 'age.end', 'status', 'subcohort') %in% colnames(ccoh.data)))>0) stop("Need the following variables in ccoh.data: ID, age.start, age.end, status, subcohort.")
  
  #######################################
  #### non-subcohort case pooling: pool based on age at diagnosis 
  #######################################
  
      ccoh.data$case <- 1*(ccoh.data$status==1 & !ccoh.data$subcohort) ## non-subcohort cases
      ccoh.case <- ccoh.data[ccoh.data$case==1,]
      ccoh.case <- ccoh.case[order(ccoh.case$age.end),]
      n.age.end <- table( ccoh.case$age.end)
      
      pool.ID <-0 ## keep track of last pool ID pools- IDs of all pools
      
      pools <- NULL #list of IDs either pooled or singles
      p.size <- NULL
      for (a in 1:length(n.age.end)){
        n.pool1 <- floor(n.age.end[a]/pool.size) # pool as much as possible
        id1<- NULL
        if (n.pool1>0) id1<- rep(1:n.pool1+pool.ID,each=pool.size)
        pool.ID <- n.pool1+pool.ID
        n.single1 <- n.age.end[a]-n.pool1*pool.size
        if (n.single1>0) id1<- c(id1,(1:n.single1)+pool.ID)
        pool.ID <- n.single1+pool.ID
        pools <- c(pools,id1)
        p.size <- c(p.size, c(rep(pool.size,n.pool1*pool.size),rep(1,n.single1)))
      }
      
      ccoh.case$pools <- pools
      ccoh.case$p.size <- p.size
      
      
      ######################################
      ##### Subcohort pooling: pool based on age at entry- cases within subcort will be followed up after dx ###
      ######################################
      ccoh.subc <-  ccoh.data[ccoh.data$case==0,]
      ccoh.subc <- ccoh.subc[order(ccoh.subc$age.start),]
      n.age.end <- table(ccoh.subc$age.start)  ## n.age.end here is actually n.age.start
      
      pools <- NULL
      p.size <- NULL
      for (a in 1:length(n.age.end)){
        n.pool1 <- floor(n.age.end[a]/pool.size) 
        id1<- NULL
        if (n.pool1>0) id1<- rep(1:n.pool1+pool.ID,each=pool.size)
        pool.ID <- n.pool1+pool.ID
        n.single1 <- n.age.end[a]-n.pool1*pool.size
        if (n.single1>0) id1<- c(id1,(1:n.single1)+pool.ID)
        pool.ID <- n.single1+pool.ID
        pools <- c(pools,id1)
        p.size <- c(p.size, c(rep(pool.size,n.pool1*pool.size),rep(1,n.single1)))
      }
      
      ccoh.subc$pools <- pools
      ccoh.subc$p.size <- p.size
      
      pool.mix <-tapply(ccoh.subc$status, ccoh.subc$pools, function(v) length(table(v))) ## check how many pools have cases in the subcohort pools
      ## pool IDs with both cases and noncases in the same pool
      pool.mixed.id <- as.numeric(names(pool.mix))[which(pool.mix>1)]
      
      ## cases in subcort age  at last fu minus one, except for the cases in a single pool
      ccoh.subc$age.end[ccoh.subc$status==1 & ccoh.subc$p.size!=1]<- ccoh.subc$age.end[ccoh.subc$status==1 & ccoh.subc$p.size!=1]-1
      status.subcohort <- rep(0,length(table(ccoh.subc$pools)))
      status.subcohort[unique(ccoh.subc$pools)%in%ccoh.subc$pools[ccoh.subc$status==1 &ccoh.subc$p.size==1]] <-1
      
      ccoh.pool <- as.data.frame(rbind(ccoh.case,ccoh.subc))
      write.table(ccoh.pool[,c('ID','pools','p.size')],row.names = F,sep='\t',quote=F,file='pool_mapping.txt')
      
      ## determine start and end age of each pool, create partial pools ##
      
      cov.mat <- NULL
      for (e1 in cov.lst) {
        cov.mat <- cbind(cov.mat,tapply(ccoh.pool[,e1], ccoh.pool$pools, sum))
      }
      rownames(cov.mat) <- names(tapply(ccoh.pool[,e1], ccoh.pool$pools, sum))
      colnames(cov.mat) <- cov.lst
      
      ### in real study pooled exposure values will be measured and saved in "Exposure_in_pools.txt" file
      ### read in exposure value and merge the exposure data with the rest of the variables
      
      age.start <- tapply(ccoh.pool$age.start,ccoh.pool$pools,max)
      age.end <- tapply(ccoh.pool$age.end,ccoh.pool$pools,min)
      
      status <- c(rep(1,length(table(ccoh.case$pools))),status.subcohort)
      subcohort <-  ccoh.pool$subcohort[!duplicated(ccoh.pool$pools)]
      pools <-  ccoh.pool$pools[!duplicated(ccoh.pool$pools)]
      p.size <-  ccoh.pool$p.size[!duplicated(ccoh.pool$pools)]
      
      ccoh.pool <- data.frame(cov.mat, age.start,age.end, status, subcohort,pools,p.size)
      
      ccoh.pool <- ccoh.pool[ccoh.pool$age.end> ccoh.pool$age.start,]
      ccoh.pool$frac <- 1
      remove( age.start,age.end, status, subcohort,pools,p.size)
      
      ####################################################################
      ## add fractional numbers for the subcohort pools with different end age due to dropout
      ####################################################################
      
      ## ccoh.subc is the subcohort individual samples
      
      ccoh.pool.org <- as.data.frame(ccoh.subc)
      n.end.age <- tapply(ccoh.pool.org$age.end,ccoh.pool.org$pools,function(v){length(table(v))})
      sum(n.end.age>1 )
      id.mixed.end <- names(n.end.age[n.end.age>1]) ## pools with different end age
      length(pool.mixed.id)
      length(id.mixed.end)
      
      partial.pool <- NULL
      
      for (id in id.mixed.end){
        pool1 <- ccoh.pool.org[ccoh.pool.org$pools==id,]
        E1 <-  as.numeric(ccoh.pool[ccoh.pool$pools==id,cov.lst])
        age.start1 <- pool1$age.start[1]
        status1 <- pool1$status[1]
        subcohort1 <- pool1$subcohort[1]
        pools1 <- pool1$pools[1]
        p.size1 <- pool1$p.size[1]
        ages <- sort(unique(pool1$age.end))
        if (length(ages)>1){
          for (j in 2:length(ages)){
            age1 <- ages[j-1]
            age2 <- ages[j]
            
            frac1 <- sum(pool1$age.end >= age2)/nrow(pool1)
            if (frac1<1 & sum(pool1$status[pool1$age.end<age2])==0) partial.pool <- rbind(partial.pool,c(E1,age1,age2,status1,subcohort1,pools1,p.size1,frac1)) 
          }
        }
      }
      
      ccoh.pool$frac <-1
      if (!is.null(partial.pool)) {
        colnames(partial.pool) <- colnames(ccoh.pool)
        ccoh.pool <- rbind(ccoh.pool, partial.pool)
      }
      ccoh.pool <- ccoh.pool[order(ccoh.pool$pools),]
      write.table(ccoh.pool,row.names = F,sep='\t',quote=F,file='pooled_data_for_analysis.txt')
    }
    
    
    ############################################################
    #### Function to call logistic regression for analysis ####
    ############################################################
    
cch.logistic <- function(ccoh.pool,cov.lst,.boot=F, n.boot=10){
      requiredPackages = c('survival','sandwich','lmtest','boot')
      for (pack1 in requiredPackages) {
        if (!require(pack1, character.only = TRUE))
          install.packages(pack1)
        library(pack1, character.only = TRUE)
      }
      
    ####################################################################
    #### create person-year format                                  ####
    ####################################################################
    
    n.rep <- ccoh.pool$age.end-ccoh.pool$age.start
    ccoh.pool2 <- ccoh.pool[rep(seq_len(nrow(ccoh.pool)), n.rep),]
    ccoh.pool2$age <-as.numeric(unlist(sapply(1:nrow(ccoh.pool),function(v){(ccoh.pool$age.start[v]+1):ccoh.pool$age.end[v]})))
    ccoh.pool2$status1 <- 0
    ccoh.pool2$status1[cumsum(n.rep)] <- ccoh.pool$status
    ### All singletons add half year to age so that they will have different stratification parameters
    ccoh.pool2$age[ccoh.pool2$p.size==1]<- ccoh.pool2$age[ccoh.pool2$p.size==1]+0.5
    ccoh.pool2<- ccoh.pool2[ccoh.pool2$subcohort|ccoh.pool2$status1==1,]
    
    #####################################################################
    ####   Logistic regression to analyze pooled case-cohort data   ####
    #####################################################################
    covs <- paste(cov.lst,collapse='+')
    
    suppressWarnings(summary(glm3.frac1<- glm(as.formula(paste("status1 ~as.factor(age)", covs,sep='+')),data=ccoh.pool2, weights=frac, family=binomial)))
    coef1 <- coeftest(glm3.frac1, vcov. = vcovCL(glm3.frac1, cluster = ccoh.pool2$pools, type = "HC0"))[cov.lst,,drop=F]
    
    
    ####################################################################
    ####    Bootstrap method to get variance                        ####
    ####################################################################
    if (.boot) {
    boot1 <- boot(unique(ccoh.pool$pools),pool.Rboot,R = n.boot, stype = "i", parallel = "multicore",ncpus=n.cpu,ccoh.pool=ccoh.pool,cov.lst=cov.lst)
    beta.boot.pool <- boot1$t
    
    est.boot1 <- NULL
    OR <- NULL
    for (k in 1:ncol(beta.boot.pool)){
    est.boot1 <- rbind(est.boot1, c(mean(beta.boot.pool[,k]),sd(beta.boot.pool[,k]),quantile(beta.boot.pool[,k],c(0.025,0.975))))

    OR=c(OR, paste(round(exp(coef1[k,1]),2)," (",round(exp(coef1[k,1]-1.96*est.boot1[2]),2),", ",round(exp(coef1[k,1]+1.96*est.boot1[2]),2),")",sep=''))
    }
    }  ## if (.boot)
    
  if (.boot) list(coef=coef1,var.est = est.boot1, OR=OR) else list(coef=coef1)
    }
    
    
    
#####################################################################################
####  Example code to call the R function for pooled sample case-cohort analysis ####  
#####################################################################################

### "ccoh_data.txt" contains the data for individual case-cohort samples 
### call create.pools to form pooled samples and data are saved in file "pooled_data_for_analysis.txt"
    read.table('ccoh_data.txt',header=T)->ccoh.data
    create.pools(ccoh.data, cov.lst=c('cov1','cov2'))

### File "Exposure_in_pools.txt" contains exposure measurements of the pooled sample. Merge with the rest of 
### the variables of the pooled samples.
    read.table('pooled_data_for_analysis.txt',header=T)->ccoh.pool
    read.table('Exposure_in_pools.txt',header=T)->E.data
    merge(ccoh.pool, E.data, by='pools')->ccoh.pool
    
### call cch.logistic function to fit the logistic regression for pooled-sample case-cohort analysis.    
    cch.logistic(ccoh.pool, c('cov1','cov2','E'),.boot=T)
