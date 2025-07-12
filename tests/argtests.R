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



