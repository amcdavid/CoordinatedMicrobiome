library(testthat)
library(mgcv)

context("Test predictions")
set.seed(123)

f = function(x, y, z){ (.5*x + .5*y)^2 + .5*z }

## 11-15-2017
# consider a sequence of sample size n = 50, 200, 2000
# p values should converge towards zero

## 11-16-2017
## generate the dataset
## n = 50    p_value = 0.4475
## n = 100   p_value = 0.5575
## n = 200   p_value = 0.5375
## n = 500   p_value = 0.895
## n = 1000  p_value = 0.975
## n = 2000  p_value = 0.995

## ??
test_that('Converge', {
  nn = c(50, 200, 800)

  p_sequence = sapply(nn, function(n){
  x = runif(n)
  y = runif(n)

  g = function(x, y){ (x^2)*y }
  Y = g(x,y) + rnorm(n, sd = 0.1)
  data_16 = as.data.frame(cbind(Y, x, y))

  ## fit the models
  fit0 = fit_si(model = Y ~ 1, x_si = cbind(x, y), gam_data = data_16)
  fit1 = gam(formula = Y ~ s(x,y), family = gaussian, data = data_16)

  ## boostrap test
  boot_test(model1 = fit1, model0 = fit0, B = 20)$p_value
  })

  expect_equal(order(p_sequence, decreasing = TRUE), 1:3)
})

## 11-17-2017
## test response name...


generate_si_data <- function(n=1000,  length_xsi = 2, Nind = 1, sigma = .1, sigmaInd = 0) {

  ## Generate mock data (Training)
  x = runif(n)
  y = runif(n)
  z = runif(n)
  u = runif(n)
  v = runif(n)

  ## Codes for individuals
  id = gl(Nind, k = 1, length = n)
  ## Individual-specific random effect
  ranef_id = rnorm(Nind, sd = sigmaInd)

  Y = f(x, y, z) + rnorm(n, sd = sigma) + ranef_id[id]
  dat = data.frame(x, y, z, u, v, Y, id)

  x_si = data.frame(x=dat$x, y=dat$y)
  if(length_xsi> 2) x_si = cbind(x_si, u=u)
  list(dat=dat, x_si=x_si, truth = f(x,y,z) + ranef_id[id], ranef_id = ranef_id)
}

train  = generate_si_data()
train3 = generate_si_data(length_xsi = 3)
test   = generate_si_data()

test_that('Can fit SI model', {
  si_fit = with(train, fit_si(model = Y ~ z, x_si = x_si, gam_data = dat))
  expect_equal(length(si_fit$alpha), 2)
  expect_equal(names(si_fit$alpha), colnames(train$x_si))
  ratio = si_fit$alpha[1]/si_fit$alpha[2]
  expect_equal(ratio, 1, tolerance = .05, check.attributes = FALSE)
  expect_is(si_fit$x_si, 'matrix')
})

test_that('Can fit SI model with random effect', {
    train_re = generate_si_data(1000, Nind = 10, sigmaInd = 1)
    si_fit = with(train_re, fit_si(model = Y ~ s(id, bs = 're'), x_si = x_si, gam_data = dat))

    ## test that fitted variance components are close to true values
})

test_that('Can predict', {
  si_fit = with(train, fit_si(model = Y ~ z, x_si = x_si, gam_data = dat))

  #fitted values/predictions are close to f(x,y,z) on training data
  si_pred = with(train, predict_si(si_fit, x_si_new = x_si, gam_data_new = dat))
  expect_equal(as.vector(si_pred), with(train$dat, f(x,y,z)), tolerance = .05,
               check.attributes = FALSE, use.names = FALSE)

  #fitted values/predictions are close to f(x,y,z) on new runif data
  si_pred1 = with(test, predict_si(si_fit, x_si_new = x_si, gam_data_new = dat))
  expect_equal(as.vector(si_pred1), with(test$dat, f(x,y,z)), tolerance = .05,
               check.attributes = FALSE, use.names = FALSE)
})

test_that('Can predict and get effects with random effect',
          {
            train_re = generate_si_data(1000, Nind = 5, sigmaInd = 5)
            si_fit = with(train_re, fit_si(model = Y ~ s(id, bs = 're'), x_si = x_si, gam_data = dat))
            si_pred = with(train_re, predict_si(si_fit, x_si_new = x_si, gam_data_new = dat))
            expect_equal(as.vector(si_pred), train_re$truth, tolerance = .05,
                         check.attributes = FALSE, use.names = FALSE)

            eff = Effect_si(1, mod = si_fit)
          })

test_that('One dimensional SI model throws an error',
  expect_error(with(train, si_fit = fit_si(model = y ~ x, x_si = as.matrix(dat$z), gam_data = dat)))
)

test_that('Can generate average effects',
{
  si_fit = with(train3, fit_si(model = Y ~ z, x_si = x_si, gam_data = dat))

  # length(focal_si) == 2
  eff = Effect_si(focal_si = c(1,2), mod = si_fit)
  expect_true(all(eff$gam_data[1,'z'] == mean(train3$dat$z)))
  predict_effect = predict_si(si_fit, eff$x_si, eff$gam_data)
  expect_equal(length(predict_effect), nrow(eff$gam_data))

  expect_equal(range(eff$x_si[,'x']), range(train3$dat$x))
  expect_equal(range(eff$x_si[,'y']), range(train3$dat$y))

# test when length(focal_si) == 1 and == 0

  # length(focal_si) == 1
  eff = Effect_si(focal_si = c(1), mod = si_fit)
  expect_true(all(eff$gam_data[1,'z'] == mean(train3$dat$z)))
  predict_effect = predict_si(si_fit, eff$x_si, eff$gam_data)
  expect_equal(length(predict_effect), nrow(eff$gam_data))

  expect_equal(range(eff$x_si[,'x']), range(train3$dat$x))
  expect_equal(range(eff$x_si[,'x']), range(train3$dat$x))

  # length(focal_si) == 0
  eff = Effect_si(focal_si = numeric(0), mod = si_fit)
  expect_true(all(eff$gam_data[1,'z'] == mean(train3$dat$z)))
  predict_effect = predict_si(si_fit, eff$x_si, eff$gam_data)
  expect_equal(length(predict_effect), nrow(eff$gam_data))

# test that focal predictors span range of training data
# added above
})
