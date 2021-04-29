library(MASS)
library(ggplot2)
library(gridExtra)
library(cowplot)

gamma=.9

testpar = c(sigma2=.02,q=0.1,A=1.5,b=.2,c0=0.1)
gen <- function(par,n){
    res = list()
    res$spikes = runif(n) < par['q']
 
    G = outer(1:n,1:n,function(x,y) gamma^(x-y))
    G[upper.tri(G)]=0

    V = gamma^(seq(0,n-1))

    res$y = par['c0']*V+par['A']*G%*%res$spikes + par['b'] + rnorm(n,0,sqrt(par['sigma2']))
    res$par=par

    return(res)
}

gibbs_sampler <- function(data,niter=100){

    # hyperparameters
    empirical_baseline.mean = mean(data[1:20])
    empirical_baseline.sd   = sd(data[1:20])
    theta.mu = c(1,empirical_baseline.mean,0)
    theta.Sigma = diag(c(1,1,1))
    theta.Sigma.inv = solve(theta.Sigma)
    alpha.sigma2=1
    beta.sigma2 = .1
    alpha.q=2
    beta.q=100

    data.len = length(data)

    par = testpar
    par['sigma2'] = 1/rgamma(n=1,shape=alpha.sigma2,rate=beta.sigma2)
    theta       = mvrnorm(n=1,mu=theta.mu, Sigma=theta.Sigma)
    while(any(theta<0)) {
        theta = mvrnorm(n=1,mu=theta.mu, Sigma=theta.Sigma)
    }

    par[c('A','b','c0')] = theta
    par['q'] = rbeta(1,alpha.q,beta.q)

    spikes = rep(0,data.len)

    par.dataframe = data.frame(as.list(par))
    spikes.mat = matrix(0,data.len,niter)

    for( i in 1:niter){

        ## flip spikes

        G = outer(1:data.len,1:data.len,function(x,y) gamma^(x-y))
        G[upper.tri(G)]=0

        V = gamma^(seq(0,data.len-1))

        W = par['A']^2*t(G)%*%G/par['sigma2']

        yt = data-par['b']-par['c0']*V

        for(t in 1:data.len){
            logalphaMH = (1-2*spikes[t])*(-sum(W[t,spikes==1]) - (1-2*spikes[t])*W[t,t]/2+par['A']/par['sigma2']*sum(G[,t]*yt)+log(par['q']/(1-par['q'])))
            if(runif(1) < exp(logalphaMH)){
                spikes[t] = 1-spikes[t]
            }
        }

        S = cbind(rowSums(G[,spikes==1]),rep(1,data.len),V)
        Lambda = solve(theta.Sigma.inv+t(S)%*%S/par['sigma2'])
        Mu = Lambda%*%(theta.Sigma.inv%*%theta.mu+t(S)%*%data/par['sigma2'])

        newpar = par
        newpar[c('A','b','c0')] = mvrnorm(n=1,mu=Mu,Sigma=Lambda)
        while(any(newpar[c('A','b','c0')]<0)){
            newpar[c('A','b','c0')] = mvrnorm(n=1,mu=Mu,Sigma=Lambda)
        }

        newpar['sigma2'] = 1/rgamma(n=1,alpha.sigma2+data.len/2,
                                        beta.sigma2+0.5*sum((data-S%*%newpar[c('A','b','c0')])^2)
                                        )
        newpar['q'] = rbeta(n=1,alpha.q+sum(spikes),beta.q+sum(1-spikes))

        par.dataframe = rbind(par.dataframe,newpar)
        spikes.mat[,i] = spikes

        par=newpar

        if(i%%10==0) cat(i,"   \r")


    }
    cat("\n")

    return(list(par=par.dataframe,spikes=spikes.mat))



} 

validation <- function(samples,gt,keep=20){
    nsamples = ncol(samples$spikes)
    prob_spikes = rowMeans(samples$spikes[,(nsamples-keep+1):nsamples])
    time = 1:length(prob_spikes)
    nspc = rowMeans(apply(samples$spikes[,(nsamples-keep+1):nsamples],2,cumsum))
    df=data.frame(time=time,
                  prob=prob_spikes,
                  nspc=nspc,
                  true_spikes=gt$spikes,
                  Y=gt$y)
    
    g_data=ggplot(df) + geom_line(aes(time,Y))+geom_vline(xintercept=time[gt$spikes==1],col="green")
    g_prob_spike = ggplot(df)+geom_line(aes(time,prob),col="red")
    g_nspc = ggplot(df)+geom_line(aes(time,nspc)) + geom_line(aes(time,cumsum(true_spikes)),col="green")

    g = plot_grid(g_data,g_prob_spike,g_nspc,ncol=1,align="v")
    plot(g)
}

gibbs.show <- function(samples,data,keep=20){
    time = 1:length(prob_spikes)
    nsamples = ncol(samples$spikes)
    prob_spikes = rowMeans(samples$spikes[,(nsamples-keep+1):nsamples])
    nspc = rowMeans(apply(samples$spikes[,(nsamples-keep+1):nsamples],2,cumsum))
    df=data.frame(time=time,
                  prob=prob_spikes,
                  nspc=nspc,
                  Y=gt$y)
    
    g_data=ggplot(df) + geom_line(aes(time,Y))
    g_prob_spike = ggplot(df)+geom_line(aes(time,prob),col="red")
    g_nspc = ggplot(df)+geom_line(aes(time,nspc)) 

    g = plot_grid(g_data,g_prob_spike,g_nspc,ncol=1,align="v")
    plot(g)
}
