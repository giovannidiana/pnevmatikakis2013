library(MASS)

gamma=.9


testpar = c(sigma2=.02,q=0.02,A=1.5,b=.2,c0=0.1)
gen <- function(par,n){
    res = list()
    res$spikes = runif(n) < par['q']
 
    G = outer(1:n,1:n,function(x,y) gamma^(x-y))
    G[upper.tri(G)]=0

    V = gamma^(seq(0,n-1))

    res$y = par['c0']*V+par['A']*G%*%res$spikes + par['b'] + rnorm(n,0,sqrt(par['sigma2']))

    return(res)
}


gibbs_sampler <- function(data,niter=100){

    # hyperparameters
    empirical_baseline.mean = mean(data[1:20])
    empirical_baseline.sd   = sd(data[1:20])
    theta.mu = c(1,empirical_baseline.mean,0)
    theta.Sigma = diag(c(1,empirical_baseline.sd^2,1))
    theta.Sigma.inv = inv(theta.Sigma)
    alpha.sigma2=1
    beta.sigma2 = .1
    alpha.q=1
    beta.q=1

    data.len = length(data)

    par = testpar
    par['sigma2'] = 1/rgamma(n=1,shape=alpha.sigma2,rate=beta.sigma2)
    theta       = mvrnorm(n=1,mu=theta.mu, Sigma=theta.Sigma)
    par[c('A','b','c0')] = theta
    par['q'] = rbeta(1,alpha.q,beta.q)

    spikes = rep(0,data.len)

    par.dataframe = data.frame(par)
    spikes.mat = matrix(0,data.len,niter)

    for( i in 1:niter){
        par.current = par

        ## flip spikes

        G = outer(1:data.len,1:data.len,function(x,y) gamma^(x-y))
        G[upper.tri(G)]=0

        V = gamma^(seq(0,data.len-1))

        W = par['A']^2*t(G)%*%G/par['sigma2']

        for(t in 1:data.len){
            logalphaMH = (1-2*spikes[t])*(-sum(W[t,spikes==1]) - (1-2*spikes[t])*W[t,t]/2+par['A']/par['sigma2']*sum(G[,t]*yt)+log(par['q']/(1-par['q'])))
            if(runif(1) < exp(logalphaMH)){
                spikes[t] = 1-spikes[t]
            }
        }

        S = cbind(par['A']*G%*%spikes,rep(1,data.len),V)
        Lambda = inv(theta.Sigma.inv+t(S)%*%S/par['sigma2'])
        Mu = Lambda%*%(theta.Sigma.inv%*%theta.mu+t(S)%*%data)

        newpar = par
        newpar[c('A','b','c0')] = mvrnorm(n=1,mu=Mu,Sigma=Lambda)
        newpar['sigma2'] = 1/rgamma(n=1,alpha.sigma2+data.len/2,
                                        beta.sigma2+0.5*sum((data-S%*%newpar[c('A','b','c0')])^2)
                                        )
        newpar['q'] = sum(spikes)/data.len

        par.dataframe = rbind(par.dataframe,newpar)
        spikes[,i] = spikes


    }



} 
