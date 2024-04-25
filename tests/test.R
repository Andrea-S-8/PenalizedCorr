library(tinytest)

set.seed(1234)
data=cbind(arima.sim(model=list(ar=0.7),n=500),arima.sim(model=list(ar=0.3),n=500))

acfdata=acf(data)
expect_equal_to_reference(acfdata,'./expected/acfdata.Rdata')

pacfdata=pacf(data)
expect_equal_to_reference(pacfdata,'./expected/pacfdata.Rdata')

ardata=ar(data)
expect_equal_to_reference(ardata,'./expected/ardata.Rdata')

# test a single lag, one series
set.seed(5923)
datama=arima.sim(model=list(ma=0.99),n=100)
acfdatama=acf(datama,lag.max=1)

# test a NND
acfdatama=acf(datama)

# test AR order larger than 3 warning
acfdatamawarn=acf(datama,estimate="invertpacf")

pacfdatama=pacf(datama)



# Colin's example
set.seed(56334)
x=arima.sim(n=1000,model=list(ar=.9))
original=acf(x,penalized=FALSE)
pen=acf(x)
invert=acf(x,estimate="invertpacf")
