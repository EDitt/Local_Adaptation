### ASTER function

# newdata= empty dataframe with variables in model
# FitnessVar = (in quotations) highest level fitness variable
# Model = model to use for predict
ASTER_predict <- function(newdata, FitnessVar, Model) {
  for (v in vars)
    newdata[[v]] <- 1
  newdata$root <- 1 #add the root
  #now reshape this dataset
  renewdata <- reshape(newdata, varying = list(vars), direction = "long", timevar = "varb", times = as.factor(vars), v.names = "resp")
  fit <- grepl(FitnessVar, as.character(renewdata$varb))
  fit <- as.numeric(fit)
  renewdata$fit <- fit
  nind <-nrow (newdata) #num rows of df
  nnode <- length (vars) 
  amat <-array(0, c(nind, nnode, nind)) #this makes an empty array
  for (i in 1:nind)
    amat[i, grep(FitnessVar, vars), i] <- 1
  predictions <- predict(Model, varvar=varb, idvar=id, root=root, newdata=renewdata, se.fit=TRUE, amat=amat)
  return(predictions)
}

