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



