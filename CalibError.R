library(riskRegression)

X1 <- qs::qread('X1')
X2 <- qs::qread('X2')

X1[["Brier"]][["score"]]
X2[["Brier"]][["score"]]

# this does not work:
# Error in KernSmooth::dpik(cumtabx/N, kernel = "box") : 
#   scale estimate is zero for input data
plotCalibration(X1,times=5,cens.method="local")
# we could not find a way to replicate this error
# and could not find where the function call KernSmooth::dpik occurs
# in PlotCalibration()...

# this works
plotCalibration(X2,times=5,cens.method="local")

# this works
plotCalibration(X1,times=5,cens.method="local",models = "SL04",bars = TRUE,
                ylim = c(0,0.05),col = c("#014667","gray"))

# this does not work:
# Error in cut.default(p, groups, include.lowest = TRUE) : 
#   'breaks' are not unique
plotCalibration(X1,times=5,cens.method="local",models = "SL05",bars = TRUE,
                ylim = c(0,0.05),col = c("#026899","gray"))

# this issue appears to arise when the 'breaks' input for the 'cut' function
# is not unique.

# in plotCalibration() the call happens with the following line:
# pcut <- cut(p, groups, include.lowest = TRUE)
# if the 'groups' created earleir inthe code is somehow not unique, we get this error

# See plotCalibraiton2.R line 338 for a quick fix; we simply set:
# groups <- unique(groups)

source('./plotCalibration2.R')

# old function (error)
plotCalibration(X1,times=5,cens.method="local",models = "SL05",bars = TRUE,
                ylim = c(0,0.05),col = c("#026899","gray"))

# new function (no error)
plotCalibration2(X1,times=5,cens.method="local",models = "SL05",bars = TRUE,
                ylim = c(0,0.05),col = c("#026899","gray"))
