library(tinytest)
# examples

##############################################
# acf.R
##############################################

### AR(1)
set.seed(1234)
data <- arima.sim(n=100, model=list(ar=0.5))

# Examples for penalized ACF/PACF and sample ACF/PACF
# penalized estimate
expect_equal_to_reference(acf(data),'./expected/exacfdata.Rdata')
# usual stats:::acf() implementation
expect_equal_to_reference(acf(data,penalized=FALSE),'./expected/exstatsacfdata.Rdata')
# penalized partial correlation estimate
expect_equal_to_reference(acf(data,type="partial"),'./expected/expacfdata.Rdata')
# usual stats:::pacf() estimate
expect_equal_to_reference(acf(data, type ="partial", penalized = FALSE),'./expected/exstatspacfdata.Rdata')
# estimate the acf by inverting the pacf
expect_equal_to_reference(acf(data,estimate="invertpacf"),'./expected/exinvertpacfdata.Rdata')

set.seed(1234)
x1 <- arima.sim(n=100, model=list(ar=0.5))
x2 <- arima.sim(n=100, model=list(ar=0.1))
x3 <- arima.sim(n=100, model=list(ar=0.9))
x <- cbind(x1, x2, x3)

#penalized estimate
expect_equal_to_reference(acf(x),'./expected/exmvacfdata.Rdata')
# usual stats:::acf() implementation
expect_equal_to_reference(acf(x, penalized = FALSE),'./expected/exmvstatsacfdata.Rdata')
# penalized partial correlation estimate
expect_equal_to_reference(acf(x, type ="partial"),'./expected/exmvpacfdata.Rdata')
# usual stats:::pacf() estimate
expect_equal_to_reference(acf(x, type ="partial", penalized = FALSE),'./expected/exmvstatspacfdata.Rdata')
# estimate the acf by inverting the pacf
expect_equal_to_reference(acf(x,estimate="invertpacf"),'./expected/exmvinvertpacfdata.Rdata')

###MA(1)
### MA(1)
set.seed(1234)
data <- arima.sim(n=100, model=list(ma=0.7))

# Examples for penalized ACF/PACF and sample ACF/PACF
# penalized estimate
expect_equal_to_reference(acf(data),'./expected/exmaacfdata.Rdata')
# usual stats:::acf() implementation
expect_equal_to_reference(acf(data,penalized=FALSE),'./expected/exmastatsacfdata.Rdata')
# penalized partial correlation estimate
expect_equal_to_reference(acf(data,type="partial"),'./expected/exmapacfdata.Rdata')
# usual stats:::pacf() estimate
expect_equal_to_reference(acf(data, type ="partial", penalized = FALSE),'./expected/exmastatspacfdata.Rdata')
# estimate the acf by inverting the pacf
expect_equal_to_reference(acf(data,estimate="invertpacf"),'./expected/exmainvertpacfdata.Rdata')

set.seed(1234)
x1 <- arima.sim(n=100, model=list(ma=0.5))
x2 <- arima.sim(n=100, model=list(ma=0.1))
x3 <- arima.sim(n=100, model=list(ma=0.9))
x <- cbind(x1, x2, x3)

#penalized estimate
expect_equal_to_reference(acf(x),'./expected/exmvmaacfdata.Rdata')
# usual stats:::acf() implementation
expect_equal_to_reference(acf(x, penalized = FALSE),'./expected/exmvmastatsacfdata.Rdata')
# penalized partial correlation estimate
expect_equal_to_reference(acf(x, type ="partial"),'./expected/exmvmapacfdata.Rdata')
# usual stats:::pacf() estimate
expect_equal_to_reference(acf(x, type ="partial", penalized = FALSE),'./expected/exmvmastatspacfdata.Rdata')
# estimate the acf by inverting the pacf
expect_equal_to_reference(acf(x,estimate="invertpacf"),'./expected/exmvmainvertpacfdata.Rdata')



##############################################
# ar.R
##############################################
set.seed(1234)
data <- arima.sim(n=100, model=list(ar=0.5))

# penalized ar model fit
expect_equal_to_reference(ar(data),'./expected/exardata.Rdata')
# default stats::ar() estimate
expect_equal_to_reference(ar(data,method="yw"),'./expected/exstatsardata.Rdata')



##############################################
# pacf.R
##############################################

### AR(1)
set.seed(1234)
data <- arima.sim(n=100, model=list(ar=0.5))

# penalized pacf
expect_equal_to_reference(pacf(data),'./expected/expacffuncdata.Rdata')
# stats::pacf
expect_equal_to_reference(pacf(data, penalized = FALSE),'./expected/exstatspacffuncdata.Rdata')


### MA(1)
set.seed(1234)
data <- arima.sim(n=100, model=list(ma=0.7))

# penalized pacf
expect_equal_to_reference(pacf(data),'./expected/exmapacffuncdata.Rdata')
# stats::pacf
expect_equal_to_reference(pacf(data, penalized = FALSE),'./expected/exmastatspacffuncdata.Rdata')

set.seed(1234)
x1 <- arima.sim(n=100, model=list(ma=0.5))
x2 <- arima.sim(n=100, model=list(ma=0.1))
x3 <- arima.sim(n=100, model=list(ma=0.9))
x <- cbind(x1, x2, x3)

# penalized pacf
expect_equal_to_reference(pacf(x),'./expected/exmvmapacffuncdata.Rdata')
# stats::pacf
expect_equal_to_reference(pacf(x, penalized = FALSE),'./expected/exmvmastatspacffuncdata.Rdata')
