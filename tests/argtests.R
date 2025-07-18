library(tinytest)


set.seed(123)
y <- arima.sim(n=100, model=list(ar=0.5))
#also checked
# y <- arima.sim(n=100, model=list(ma=0.5))
y1 <- c("a", "b", "c", "d", "e", "f")
y2 <- c("a", "1", "c", "3", "e", "6")


#testing that PenalizedCorr can only be logical
### acf ###
expect_error(acf(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(acf(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(acf(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(acf(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### pacf ###
expect_error(pacf(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(pacf(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(pacf(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(pacf(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### ar ###
expect_error(ar(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(ar(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### auto.acf ###
expect_error(auto.acf(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(auto.acf(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(auto.acf(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(auto.acf(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### invertpacf ###
expect_error(invertpacf(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(invertpacf(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(invertpacf(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(invertpacf(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### DLpencoaf ###
expect_error(DLpencoef(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(DLpencoef(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(DLpencoef(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(DLpencoef(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected


### ar.penyw ###
expect_error(ar.penyw(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar.penyw(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar.penyw(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(ar.penyw(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected






##testing that data series can only be numeric
#### acf ####
expect_error(acf(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(acf(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(acf(y2), "'x' must be numeric")
#Passed as error found as expected

#### pacf ####
expect_error(pacf(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(pacf(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(pacf(y2), "'x' must be numeric")
#Passed as error found as expected

#### auto.acf ####
expect_error(auto.acf(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(auto.acf(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(auto.acf(y2), "'x' must be numeric")
#Passed as error found as expected

#### invertpacf ####
expect_error(invertpacf(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(invertpacf(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(invertpacf(y2), "'x' must be numeric")
#Passed as error found as expected

#### DLpencoef ####
expect_error(DLpencoef(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(DLpencoef(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(DLpencoef(y2), "'x' must be numeric")
#Passed as error found as expected

#### ar ####
expect_error(ar(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(ar(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(ar(y2), "'x' must be numeric")
#Passed as error found as expected

#### ar.penyw ####
expect_error(ar.penyw(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(ar.penyw(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(ar.penyw(y2), "'x' must be numeric")
#Passed as error found as expected





set.seed(12345)
x1 = arima.sim(model=list(ar=0.4),n=500)
x2 = arima.sim(model=list(ar=0.8),n=500)
x=cbind(x1, x2)

#testing that PenalizedCorr can only be logical
### acf ###
expect_error(acf(x,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(acf(x,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(acf(x,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(acf(x,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### pacf ###
expect_error(pacf(x,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(pacf(x,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(pacf(x,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(pacf(x,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### ar ###
########
#ar doesn't find penalized so all failed even when should have passed
########
expect_error(ar(x,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar(x,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar(x,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(ar(x,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### auto.acf ###
########
#auto.acf doesn't find penalized so first is failed when should have passed
########
expect_error(auto.acf(x,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(auto.acf(x,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(auto.acf(x,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(auto.acf(x,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### invertpacf ###
########
#invertpacf doesn't find penalized so first is failed when should have passed
########
expect_error(invertpacf(x,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(invertpacf(x,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(invertpacf(x,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(invertpacf(x,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### DLpencoaf ###
expect_error(DLpencoef(x,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(DLpencoef(x,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(DLpencoef(x,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(DLpencoef(x,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### ar.penyw ###
########
#ar doesn't find penalized so all failed even when should have passed
########
expect_error(ar.penyw(x,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar.penyw(x,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar.penyw(x,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(ar.penyw(x,penalized = 3), "penalized must be logical")
#Passed as errors found as expected




##testing that data series can only be numeric
#### acf ####
################
#fails but not for the right reasons
################
expect_error(acf(x), "'x' must be numeric")
#Failed as no error found as expected

#### pacf ####
expect_error(pacf(x), "'x' must be numeric")
#Failed as no error found as expected

#### auto.acf ####
################
#fails but not for the right reasons
################
expect_error(auto.acf(x), "'x' must be numeric")
#Failed as no error found as expected

#### invertpacf ####
################
#fails but not for the right reasons
################
expect_error(invertpacf(x), "'x' must be numeric")
#Failed as no error found as expected

#### DLpencoef ####
expect_error(DLpencoef(x), "'x' must be numeric")
#Failed as no error found as expected

#### ar ####
################
#fails but not for the right reasons
################
expect_error(ar(x), "'x' must be numeric")
#Failed as no error found as expected

#### ar.penyw ####
################
#fails but not for the right reasons
################
expect_error(ar.penyw(x), "'x' must be numeric")
#Failed as no error found as expected




set.seed(12345)
m <- arima.sim(n=100, model=list(ar=0.8))
#expect_equal(m, readRDS("./tests/expected/equals-error/m"))


