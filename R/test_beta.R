apow <- function(x,n){
    if(n==0){
        return(1)
    } else {
        return(x*apow(x+1,n-1))
    }
}

marg_prob <- function(n,t,a,b){
    return(apow(a,n)*apow(b,t-n)/apow(a+b,t))
}
