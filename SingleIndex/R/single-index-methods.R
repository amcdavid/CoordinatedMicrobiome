# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

si <- function(theta, model, x_si, gam_data, opt = TRUE, fx = FALSE, method = 'ML') {
  ## Fit single index model using gam call, given theta (defines alpha).
  ## Return ML if opt==TRUE and fitted gam with theta added otherwise.
  ## Suitable for calling from 'optim' to find optimal theta/alpha.
  alpha <- c(1,theta) ## constrained alpha defined using free theta
  names(alpha) = colnames(x_si)
  x_si = as.matrix(x_si)
  stopifnot(length(alpha) == ncol(x_si))
  kk <- sqrt(sum(alpha^2))
  alpha <- alpha/kk  ## so now ||alpha||=1
  a <- x_si %*% alpha     ## argument of smooth
  gam_data = cbind(gam_data, fx = fx, a = a)
  model1 = update.formula(model, . ~ . + s(a, fx = fx, k = 20))
  # so we can find variables defined here
  environment(model1) = environment()

  b <- gam(model1, family=gaussian, method=method, data = gam_data) ## fit model
  if (opt) return(b$gcv.ubre) else {
    b$alpha <- alpha  ## add alpha
    J <- outer(alpha, -theta/kk^2) ## compute Jacobian
    for (j in 1:length(theta)) J[j+1, j] <- J[j+1, j] + 1/kk
    b$J <- J ## dalpha_i/dtheta_j
    return(b)
  }
}

fit_si = function(model, x_si, gam_data){
  #fit the single index using the three lines of code below
  #return an object of class gam
  if(ncol(x_si)<2) stop('Not supported for univariate single index models; call GAM directly.')
  init_theta = rep(1, ncol(x_si)-1)
  if(is.null(colnames(x_si))) colnames(x_si) = paste0('si', '.', seq_len(ncol(x_si)))
  f0 <- nlm(si, init_theta, model, x_si, gam_data, fx = TRUE)
  f1 <- nlm(si, f0$estimate, model, x_si, gam_data, hessian = TRUE, fx = FALSE)
  b <- si(f1$estimate, model, x_si, gam_data, opt = FALSE, method = 'GCV.Cp')
  b$x_si = as.matrix(x_si)
  b$gam_data = gam_data
  return(b)
}

# ## function boot_si will generate one bootstrap sample about the single index model
boot_si = function(si_fit) {
 # r = b$gam_data$y - predict(b)
 r = residuals(si_fit, type = 'response')
 value = fitted(si_fit) + sample(x = r, size = length(r), replace = T)
 gam_data_boot = si_fit$gam_data
 response_name = as.character(formula(si_fit))[[2]]
 gam_data_boot[[response_name]] = value
 return(gam_data_boot)
}

boot_test = function(model1, model0, B, return_model = FALSE) {
  # model1 for GAM
  # model0 for fit_si
  # B is the number of resampling
  this_ll = rep(NA, times = B)
  boot_model = list()

  for (b in seq_len(B)) {
    boot_data = boot_si(model0)
    m_boot = update(model1, data = boot_data)
    if(return_model) boot_model[[b]] = m_boot
    this_ll[b] = logLik(m_boot)
  }

  obs_ll = logLik(model1)
  p_value = mean(this_ll > obs_ll)
  return(list(p_value = p_value, boot_ll = this_ll, model1_ll = obs_ll, boot_model = boot_model))
}

predict_si = function(gam_si, x_si_new, gam_data_new, ...){

  # check that x_si_new has same number of rows as gam_data_new
  if (nrow(x_si_new) != nrow(gam_data_new))
    stop("The x_si_new should have the same number of rows as the gam_data_new.")

  # coerce x_si_new to matrix
  x_si_new = as.matrix(x_si_new)

  # check that x_si_new has the right column names/order (matches alpha)
  if (all(colnames(x_si_new) != names(gam_si$alpha)))
    stop("The x_si_new should have the same column names/order as the alpha.")

  a_new = x_si_new %*% gam_si$alpha
  gam_data_new$a = a_new
  predict(gam_si, newdata = gam_data_new, ...)

}

# for testing
# test  = generate_si_data(length_xsi = 3)
# mod = with(test, fit_si(model = Y ~ z, x_si = x_si, gam_data = dat))
# focal_si = 1

typical_value = function(matrix){
  # return colMean if numeric
  # return mode if factor
  n_col = ncol(matrix)
  typical_v = list()

  for (i in seq_len(n_col)) {
    if(is.factor(matrix[,i])){
      col = matrix[,i]
      tab = table(col)
      typical_v[[i]] = factor(levels(col)[which.max(tab)], levels(col))
    }
    else typical_v[[i]] = mean(matrix[,i])
  }
  names(typical_v) = colnames(matrix)
  return(typical_v)
}

# dataframe = expand.grid(c('A', 'B'), 1:5)
# constant_list = list(one=1, two=2, three=factor('A'))
# cbind(dataframe, constant_list)

Effect_si = function(focal_si, focal_gam = numeric(0), mod,
                     length_out = rep(20, length = length(focal_si)), ...){
  # focal_si gives column indices into mod$x_si
  # focal_gam is unsupported for the moment

  if(length(focal_gam) > 0) stop("Unsupported")
  stopifnot(length(length_out)==length(focal_si))


  ## browser(expr = length(focal_si)==0)
  focal_si = seq_len(ncol(mod$x_si)) %in% focal_si
  focal_gam = seq_len(ncol(mod$model)) %in% focal_gam

  # for focal_si, get min and max of relevant columns in mod$x_si
  focal_si_min = apply(mod$x_si[,focal_si, drop = FALSE], 2, min)
  focal_si_max = apply(mod$x_si[,focal_si, drop = FALSE], 2, max)

  # for everything else (fixed_predictors), take the average of the design matrix
  other_data = mod$x_si[,!focal_si, drop = FALSE]
  other_mean = typical_value(other_data) ## ??

  # and vary the focal predictors over their ranges
  x_si_pred_focal_list = lapply(seq_along(focal_si_min), function(i) {
    seq(focal_si_min[i], focal_si_max[i], length.out = length_out[i])
  })
  names(x_si_pred_focal_list) = colnames(mod$x_si)[focal_si]
  x_si_pred_focal = do.call(expand.grid, x_si_pred_focal_list)

  # generate x_pred, x_si_pred which hold the fixed_predictors at their fixed values
  if(nrow(x_si_pred_focal)>0 & length(other_mean)>0){
    ## repeat the elements of the list as the dataframe
    x_si_pred = cbind(x_si_pred_focal, other_mean)
  } else if (length(other_mean)>0) {
    x_si_pred = as.data.frame(other_mean)
  } else if (nrow(x_si_pred_focal)>0) {
    x_si_pred = as.data.frame(x_si_pred_focal)
  } else{
    stop('Illegal path.')
  }

  # return x_pred, x_si_pred and predicted values
  other_data_x = mod$model[,!focal_gam, drop = FALSE]
  other_mean_x = typical_value(other_data_x) ##
  all_pred = cbind(x_si_pred, other_mean_x)
  pred = predict_si(mod, x_si =  all_pred[,colnames(mod$x_si)],
                    gam_data = all_pred[,colnames(mod$model)], se.fit = TRUE)
  cbind(.fitted = pred, all_pred)
}
