compute.sampling.multiplier=function(model,data,lower=0,upper=90)
{
   TotalArea=integrate(integrate.gam,lower=lower,upper=upper,model=model,vis=1,
            beauf=0)$val
   int.areas=apply(data[,c("begin","end","vis","beaufort")],1,function(x)
            integrate(integrate.gam,lower=x[1],upper=x[2],model=model,vis=as.vector(x[3]),
            beauf=as.vector(x[4]))$val)
   return(TotalArea/sum(int.areas))
}
integrate.gam=function(x,model,vis,beauf)
{
# Function for computing area under portions of the gam migration model
#
# Arguments:
#    x     -  vector of times to be evaluated within range of migration curve
#  model   -  gam model fitted to migration curve data
#  vis     - vis values to use for prediction if needed
#  beauf   - beaufort values to use for prediction if needed
#
 newdata=data.frame(time=x,vis=vis,beauf=beauf)
 return(predict(model,newdata=newdata,type="response"))
}
