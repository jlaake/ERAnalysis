io.glm <-
function(datavec, fitformula, eps = 0.00001, iterlimit = 500)
{
# ---------------------------------------------------------------
#  This is the code that uses the iterative offset glm or gam
#  approach; iteration is done until parameters are within a
#  certain epsilon (eps) or iteration limit (iterlimit) exceeded.
#
#  Note: David used offset in formula and I've put it as an
#  argument to glm and gam functions.
#
# Input : datavec = dataframe
#         fitformula = formula
#         eps = convergence criterion - fixed
#         iterlimit = maximum number of iterations allowed - fixed
#
# Output: list with
#       Note: modified to return glm object only with
#               class("ioglm","glm","lm")
#               class("ioglm","gam")
#
# $glmobj:  glm model
# $offsetvalue: final offsetvalues from iterative fit
# $plotobj: gam plot object (if GAM & gamplot==TRUE, else NULL)
# ----------------------------------------------------------------
#
  done <- FALSE
  i <- 1
  plotobj <- NULL
  if(is.null(datavec$offsetvalue))datavec$offsetvalue=0
  while(i <= iterlimit & !done) {
#  fit the glm or gam
    offsetvalue=datavec$offsetvalue
    ioglm <- glm(formula = fitformula, family = binomial, data = datavec, offset=offsetvalue)
    coeff <- ioglm$coeff
    fittedp <- ioglm$fitted.values

    if(i == 1) {
      oldmodel <- ioglm
      oldcoeff <- coeff
      oldp <- fittedp
    }else{
#    calculate differences between previous and present set of model outputs
      reldiff <- max(abs(plogis(coeff) - plogis(oldcoeff))/plogis(oldcoeff))

      if(is.na(reldiff)) {
        print("Can't calculate regression coefficients - model has not converged")
        print(" - last fit used for estimation" )
        ioglm <- oldmodel
        done <- TRUE
      }

      if(reldiff < eps & !done) {
        done <- TRUE
      }else{
        oldmodel <- ioglm
        oldcoeff <- coeff
        oldp <- fittedp
      }
    }
    if(!done){
      oldoff <- datavec$offsetvalue
        off <-  - log(plogis(predict(ioglm) - datavec$offsetvalue))
      datavec$offsetvalue[datavec$station == "S"] <- off[
        datavec$station == "P"]
      datavec$offsetvalue[datavec$station == "P"] <- off[
        datavec$station == "S"]
    }
    i <- i + 1
  }

  if(!done){
    datavec$offsetvalue <- oldoff
    warning("Iteration limit exceeded - last fit used for estimation")
  }

  class(ioglm)=c("ioglm",class(ioglm))

  return(ioglm)
}

