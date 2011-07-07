fit.fS <-
function(par,pstable,gsS,nmax=20)
{
# Function to fit observed pod size distribution and estimate true
# pod size distribution based on gamma distribution

  ps=gammad(par,nmax)
  -sum(pstable*log(colSums(gsS*ps)))
}

