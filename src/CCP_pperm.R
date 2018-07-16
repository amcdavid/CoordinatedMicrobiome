fit_cancor = function(X, Y, xcenter, ycenter){
    cc = cancor(X, Y, xcenter = xcenter, ycenter = ycenter)
    cc$Xdat = X
    cc$Ydat = Y
    cc
}

# fit_cancor = function(X, Y, xcenter, ycenter){
#     cc = PMA::CCA(X, Y, xcenter = xcenter, ycenter = ycenter)
#     cc$Xdat = X
#     cc$Ydat = Y
#     cc
# }

eval_cancor = function(fit, Xnew, Ynew){
    stopifnot(inherits(fit, 'FittedCancor'))
    if(missing(Xnew)){
        Xnew = fit$Xdat
    }
    if(missing(Ynew)){
        Ynew = fit$Ydat   
    }
    if(!is.null(fit$xcenter)){
        Xnew = t(t(Xnew - fit$xcenter))
    }
    if(!is.null(fit$ycenter)){
        Ynew = t(t(Ynew - fit$ycenter))
    }
    xscore = Xnew %*% fit$xcoef
    yscore = Ynew %*% fit$ycoef
    p1 = ncol(Xnew)
    p2 = ncol(Ynew)
    min_p = min(p1, p2)
    cxy = cor(xscore, yscore)[1:min_p, 1:min_p, drop = FALSE]
    energy = diag(cxy)^2 - (sum(cxy^2) - sum(diag(cxy)^2))/(min_p*(min_p-1))
    list(cor = diag(cxy), energy = cumsum(diag(cxy)^2)/sum(cxy^2))
}

fit_lm = function(Y, data_, formula_){
    data_ = as.data.frame(data_)
    fx  = as.formula(paste(c('Y', as.character(formula_)), collapse = ''))
    lmx = lm(fx, data=data_)
    lmx$Y = Y
    lmx$data_ = data_
    lmx
}

resid_lm = function(fit, Y, data_){
    if(missing(Y)) Y = fit$Y
    if(missing(data_)) data_ = fit$data_
    data_ = as.data.frame(data_)
    Yhat = predict(fit, newdata = data_)
    Y - Yhat
}

do_cc = function(xbiom, ybiom, formula_, center){
    fx = fit_lm(getbiom(xbiom), colData(xbiom), formula_)
    rx = resid(fx)
    fy = fit_lm(getbiom(ybiom), colData(xbiom), formula_)
    ry = resid(fy)
    fcc = fit_cancor(rx, ry, xcenter = center, ycenter = center)
    fcc$formula_ = formula_
    fcc$fx = fx
    fcc$fy = fy
    class(fcc) = c('FittedCancor', class(fcc))
    fcc
}

crossval_cc = function(xbiom, ybiom, nfolds = 10, ...){
    xdat = colData(xbiom) %>% as.data.frame() %>% mutate(idx = seq_len(nrow(.)))
    strata = xdat$Participant.ID %>% unique()
    folds = data_frame(Participant.ID = sample(strata), fold = rep(seq_len(nfolds), length.out = length(strata)))
    
    cv_folds = purrr::map(seq_len(nfolds), function(i){
        trainidx = right_join(xdat, folds %>% filter(fold != i), by = 'Participant.ID')$idx
        dcc = do_cc(xbiom[,trainidx], ybiom[,trainidx], ...)
        e_train = eval_cancor(dcc)
        e_test = predict_FittedCancor(dcc, xbiom[,-trainidx], ybiom[,-trainidx])
        list(train = e_train, test = e_test)
    })
    
    ## Pull out CV metrics
    outputs = expand.grid(sample = c('test', 'train'), metric = c('cor', 'energy'), stringsAsFactors = FALSE)
    outputs = outputs %>% group_by_all() %>% do( fold_mat = map(cv_folds, pluck, .$sample, .$metric) %>% do.call('rbind', .)) 
    ## %>% mutate(mean = list(colMeans(matrix)), sd = rowSds(matrix))
    outputs = outputs %>%ungroup() %>% mutate(mean = map(.$fold_mat, ~colMeans(.x)), 
                                              sd = map(.$fold_mat, ~colSds(.x)),
                                              coordinate = map(.$fold_mat, ~seq_len(ncol(.x)))) %>% select(-fold_mat)
    outputs %>% unnest(mean, sd, coordinate)
}


predict_FittedCancor = function(object, xbiom, ybiom, ...){
    if(missing(xbiom)){
        X = object$fcc$Xdat
        X_data = object$lmx$data_
    } else{
        X = getbiom(xbiom)
        X_data = colData(xbiom)
    }
    if(missing(ybiom)){
        Y = object$fcc$Ydat
        Y_data = object$lmy$data_
    } else{
        Y = getbiom(ybiom)
        Y_data = colData(ybiom)
    }
    
    X = resid_lm(object$fx, X, X_data)
    Y = resid_lm(object$fy, Y, Y_data)
    
    eval_cancor(object, X, Y)
    
}

get_resids = function(list_of_se, formula_){
    x = getbiom(list_of_se[[1]])
    y = getbiom(list_of_se[[2]])
    fx  = as.formula(paste(c('x', as.character(formula_)), collapse = ''))
    fy  = as.formula(paste(c('y', as.character(formula_)), collapse = ''))
    lmx = lm(fx, data=colData(list_of_se[[1]]))
    lmy = lm(fy, data=colData(list_of_se[[2]]))
    x = resid(lmx)
    y = resid(lmy)
    Hmisc::llist(x,y, lmx, lmy)
}


## Obsolete code follows

# From CCP::pperm
# Updated to allow parallel evalution of rhostart
ccp_p.perm = function (X, Y, nboot = 999, rhostart = seq_len(min(ncol(X), ncol(Y))), type = "Wilks") {
    if (nrow(Y) != nrow(X)) 
        stop(" Function p.perm: Input arrays must not have unequal number of rows.\n")
    if ((type == "Roy") && (rhostart > 1)) 
        stop(" Function p.perm: When using Roy's Largest Root, parameter rhostart can only be 1.\n")
    N <- nrow(Y)
    ind <- 1:N
    p = dim(X)[2]
    q = dim(Y)[2]
    minpq = min(p, q)
    if (any(rhostart > minpq))
        stop(" Function p.perm: Parameter rhostart too big.\n")
    rho0 <- cancor(X, Y[ind, ])$cor
    if (type == "Wilks") 
        stat0 <- WilksLambda(rho0, p, q)[rhostart]
    else if (type == "Hotelling") 
        stat0 <- HotellingLawleyTrace(rho0, p, q)[rhostart]
    else if (type == "Pillai") 
        stat0 <- PillaiBartlettTrace(rho0, p, q)[rhostart]
    else if (type == "Roy") 
        stat0 <- p.Roy(rho0, N, p, q)$stat
    else stop(" Function p.perm: Illegal type of test statistic.\n")
    stat <- matrix(NA, nrow = length(rhostart), ncol = nboot)
    for (i in 1:nboot) {
        ind.mixed <- sample(ind, size = N, replace = FALSE)
        rho <- cancor(X, Y[ind.mixed, ])$cor
        if (type == "Wilks") 
            stat[,i] <- WilksLambda(rho, p, q)[rhostart]
        else if (type == "Hotelling") 
            stat[,i] <- HotellingLawleyTrace(rho, p, q)[rhostart]
        else if (type == "Pillai") 
            stat[,i] <- PillaiBartlettTrace(rho, p, q)[rhostart]
        else if (type == "Roy") 
            stat[,i] <- p.Roy(rho, N, p, q)$stat
    }
    if (type == "Wilks") {
        nexcess <- rowSums(stat <= stat0)
        p <- rowMeans(stat <= stat0)
    }
    else {
        nexcess <- rowSums(stat >= stat0)
        p <- rowMeans(stat >= stat0)
    }
    invisible(list(id = "Permutation", type = type, stat0 = stat0, 
                   stat = stat, nexcess = nexcess, p.value = p))
}

# mlr code for canonical correlation
makeRLearner.regr.cancor = function(){
    makeRLearnerRegr('regr.cancor',
                     package='stats',
                     properties = 'numerics',
                     par.set = makeParamSet(makeIntegerLearnerParam("dim_x", lower=1),
                                            makeIntegerLearnerParam("dim_rho", lower=1)),
                     name = 'Canonical Correlation Analysis',
                     short.name = 'cancor')
}

#We have to glue together X and Y matrices.
#Pull them apart using the dim_x parameter
getXY = function(d, dim_x){
    if(dim_x < 1 || dim_x > ncol(d)) stop('dim_x out of bounds')
    ix = seq_len(dim_x)
    list(x=as.matrix(d[,ix]),
         y=as.matrix(d[,-ix]))
}

#Fit canonical correlation 
trainLearner.regr.cancor = function(.learner, .task, .subset, .weights = NULL, dim_x, dim_rho, ...){
    d = mlr::getTaskData(.task, .subset)
    xy = getXY(d, dim_x)
    x = xy$x
    y = xy$y
    cc = cancor(x, y)
    cc$min_pq = min(ncol(x), ncol(y))
    cc
}

# Ignores type and component arguments currently
predictLearner.regr.cancor = function(.learner, .model, .newdata, type='rho', component=1L, ...) {
    xc = .model$learner.model$xcoef
    yc = .model$learner.model$ycoef
    min_pq = min(.model$learner.model$min_pq, nrow(.newdata))
    dim_x = .model$learner$par.vals$dim_x
    xy = getXY(.newdata, dim_x)
    xscore = xy$x %*% xc
    yscore = xy$y %*% yc
    if (FALSE) {
        type = match.arg(type, c('x', 'y', 'rho'))
        if (length(component) != 1 || !is.numeric(component)) stop('`component` must be scalar')
        if (component > min_pq) stop('`component` exceeds minimum number of columns in `X` and `Y`')
        dlist = list(...)
        if (type == 'x') return(xscore[,component])
        if (type == 'y') return(yscore[,component])
    }
    out = rep(NA, nrow(xy$x))
    out[seq_len(min_pq)] = diag(cor(xscore, yscore))[seq_len(min_pq)]
    return(out)
}

cancor.rho.fun = function(task, model, pred, feats, extra.args){
    components = 1
    if (!is.null(extra.args$components)) components = extra.args$components
    min(pred$data[components,2])
}

test.mean = mlr::test.mean
cancor.rho1 = mlr::makeMeasure('cancor.rho', minimize = FALSE, properties = c('regr', 'req.model'), fun = cancor.rho.fun,
            extra.args = list(components=1), aggr = test.mean, best = 1, worst = -1)


cancor.rho2 = mlr::makeMeasure('cancor.rho', minimize = FALSE, properties = c('regr', 'req.model'), fun = cancor.rho.fun,
                               extra.args = list(components=1:2), aggr = test.mean, best = 1, worst = -1)

cancor.rho3 = mlr::makeMeasure('cancor.rho', minimize = FALSE, properties = c('regr', 'req.model'), fun = cancor.rho.fun,
                               extra.args = list(components=1:3), aggr = test.mean, best = 1, worst = -1)

cancor.rho4 = mlr::makeMeasure('cancor.rho', minimize = FALSE, properties = c('regr', 'req.model'), fun = cancor.rho.fun,
                               extra.args = list(components=1:4), aggr = test.mean, best = 1, worst = -1)



