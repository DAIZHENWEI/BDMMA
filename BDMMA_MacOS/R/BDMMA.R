#' Bayesian Dirichlet--Multinomial approach for meta-analysis of metagenomic read counts
#' @import Rcpp RcppArmadillo RcppEigen
#' @param X A data.frame for covariates, including main effect variable and confounding variables.
#' The first column represents the main effect variable.
#' @param Y A data.frame for taxonomic read counts with the rows representing clades and columns representing taxa.
#' @param batch A numeric vector labeling the batch for each sample.
#' @param abundance_threshold The minimum abundance level for the taxa to be included (default value = 5e-05).
#' @param burn_in The length of burn in period before sampling the parameters (default value = 5,000).
#' @param sample_period The length of sampling period for estimating parameters' distribution (default value = 5,000)
#' @param bFDR The false discovery rate level to control (default value = 0.1).
#' @return A list contains the selected taxa and summary of parameters included in the model.
#' \item{selected.taxa}{A list includes the selected taxa fesatures that are significantly associated
#' with the main effect variable.}
#' \item{parameter_summary}{A data.frame contains the mean and quantiles of the parameters included
#' in the model. Each row includes a parameter's distribution summary and the parameter name is
#' labeled in the first row. alpha_g: the baseline intercept of g-th taxon; betaj_g: the association strength between
#' the g-th taxon and j-th input variables; deltai_g: the batch effect parameter of batch i, taxon g;
#' L_g: the posterior selection probability of g-th taxon; p: the proportion of significantly associated
#' taxa; eta: the standard deviation of the spike distribution (in the spike-and-slab prior).}
#' @docType data
#' @examples
#' data(dat)
#' X <- dat$X
#' Y <- dat$Y
#' batch <- dat$batch
#' continuous <- dat$continuous
#' output <- BDMMA(X, Y, batch, continuous, burn_in = 3000, sample_period = 3000)
#' @export

BDMMA=function(X, Y, batch, continuous, abundance_threshold = 0.00005, burn_in = 5000,
               sample_period = 5000, bFDR = 0.1){

  batch = as.numeric(factor(batch))

  # check input data
  if (length(unique(X[,1]))>2|length(unique(X[,1]))<2){
    stop("The main effect variable is not binary")
  }

  if (length(unique(batch))<2){
    stop("Batch number is less than two: not a well defined batch effect indicator")
  }

  # filter low abundance bacteria
  abundance_m = sweep(Y, 1, rowSums(Y), "/")
  taxa_abundance = apply(abundance_m, 2, mean)
  YY = Y[,(taxa_abundance>abundance_threshold)]
  s_prop = taxa_abundance[taxa_abundance>abundance_threshold]
  if (length(s_prop)<ncol(Y)){
    YY = cbind(YY, (rowSums(Y)-rowSums(YY)))
    s_prop = c(s_prop,(1-sum(s_prop)))
  }


  # data characters
  n_dim = ncol(YY)
  n_size = rowSums(Y)
  n_var = ncol(X)
  n_batch = length(unique(batch))
  taxa = names(YY)

  # Normalize clinical metadata
  X = sweep(X, 2, colMeans(X), "-")
  sd_X=apply(X, 2, sd)
  if (sum(continuous != 0) > 0){
    X[,(continuous != 0)] = as.matrix(X[,(continuous != 0)])
    X[,(continuous != 0)] = sweep(X[,(continuous != 0)], 2, sd_X[(continuous != 0)], "/")
  }

  # Initialize the parameter
  a = b = 1
  c = 0.1
  sigma1 = sigma2 = sigma3 = sqrt(10)
  eta = 0.1
  L = alpha = rep(0, n_dim)

  lambda = rexp(c, n = 1)
  p = rbeta(1, shape1 = 1, shape2 = 1)

  for (i in 1:n_dim){
    alpha[i] = rexp(1/(lambda*s_prop[i]), n = 1)
    L[i] = rbinom(size = 1, n = 1, prob = p)
  }

  alpha_m = matrix(alpha, nrow = nrow(YY), ncol = n_dim, byrow = T)
  beta = matrix(0, nrow = n_var, ncol = n_dim)
  delta = matrix(0, nrow = n_batch, ncol = n_dim)
  delta_m = matrix(0, nrow = nrow(YY), ncol = n_dim)

  # Estimate parameters with MCMC
  iter = burn_in + sample_period

  X = as.matrix(X)
  YY = as.matrix(YY)

  cat("#################### Start MCMC ####################\n\n")

  output = .Call('_BDMMA_MCMC', PACKAGE = 'BDMMA', alpha = alpha, alpha_m = alpha_m, x = X, y = YY,
                beta = beta, delta = delta, delta_m = delta_m, e_delta = t(delta),
                T = n_dim, N = nrow(YY), K = n_var, I = n_batch, lambda = lambda,
                prop = s_prop, L = L, sigma1 = sigma1, sigma2 = sigma2, sigma3 = sigma3,
                iter = iter, eta = eta, a = a, b = b, p = p, batch = batch)

  ## calculate the posterior mean of parameters
  alpha_mean = apply(output[(iter-sample_period):iter, 1:n_dim], 2, mean)
  beta1_mean = apply(output[(iter-sample_period):iter, (n_dim+1):(2*n_dim)], 2, mean)
  L_mean = apply(output[(iter-sample_period):iter, (2*n_dim+1):(3*n_dim)], 2, mean)
  eta_mean = mean(output[(iter-sample_period):iter, (3*n_dim+1)])
  p_mean = mean(output[(iter-sample_period):iter, (3*n_dim+2)])

  if (n_var >= 2){
    beta2_mean = apply(output[(iter-sample_period):iter, (3 * n_dim + 3):((2 + n_var) * n_dim + 2)], 2, mean)
  }else{beta2_mean = NULL}

  delta_mean = apply(output[(iter-sample_period):iter, ((2 + n_var) * n_dim + 3):((2 + n_var + n_batch) * n_dim + 2)], 2, mean)

  quantile_2 = function(x){
    return(quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975)))
  }

  ## Generate the summary table for the parameters
  parameter_summary = t(apply(output[(iter-sample_period):iter,], 2, quantile_2))
  parameter_summary = data.frame(mean = c(alpha_mean, beta1_mean, L_mean, eta_mean, p_mean,
                                          beta2_mean, delta_mean), parameter_summary)

  name1 = name2 = name3 = rep(0, n_dim)
  name4 = matrix(0, nrow = n_var - 1, ncol = n_dim)
  name5 = matrix(0, nrow = n_batch, ncol = n_dim)

  for (ii in 1:n_dim){
    name1[ii] = paste("alpha", ii, sep = "_")
    name2[ii] = paste("beta1", ii, sep = "_")
    name3[ii] = paste("L", ii, sep = "_")
    if (n_var>=2){
      for (jj in 2:n_var){
        name4[jj-1,ii] = paste(paste("beta", jj, sep=""), ii, sep = "_")
      }
    }
    else name4 = NULL
    for (kk in 1:n_batch){
      name5[kk,ii] = paste(paste("delta", kk, sep=""), ii, sep = "_")
    }
  }
  name4 = as.vector(name4)
  name5 = as.vector(name5)
  row.names(parameter_summary) = c(name1, name2, name3, "eta", "p", name4, name5)

  ## Record the trace of the parameters
  trace = as.data.frame(output)
  names(trace) = c(name1, name2, name3, "eta", "p", name4, name5)

  ## Select the significantly associated taxa
  prediction_1 = (L_mean > 0.5) * 1
  cutoff = fdr_cut(L_mean, alpha = bFDR)
  prediction_2 = (L_mean >= cutoff) * 1
  selected.taxa = list()
  selected.taxa$MIM = taxa[prediction_1 > 0]
  selected.taxa$bFDR = taxa[prediction_2 > 0]

  output=list()

  output$trace = trace
  output$parameter_summary = parameter_summary
  output$selected.taxa = selected.taxa

  return(output)
}


#' Visualize batch effect with principal coordinate analysis
#' @param Y Y A data.frame for taxonomic read counts with the rows representing clades
#' and columns representing taxa.
#' @param batch A numeric vector labeling the batch for each sample.
#' @param main_variable Optional. A vector containing the main effect variable. Only for
#' categorical main effect variable.The function will generate a figure for each catagory.
#' @param method A string indicating which method should be used to calculate the distance matrix for
#' principal coordinate analysis.
#' @return The function returns a list containing plot objects of principal coordinate analysis figures.
#' @examples
#' data(dat)
#' figure <- VBatch(dat$Y, batch = dat$batch, main_variable = dat$X[,1], method = "bray")
#' print(figure[[1]])
#' print(figure[[2]])
#' @export
VBatch = function(Y, batch, main_variable = NULL, method = "bray"){
  batch = as.factor(batch)
  nsize = rowSums(Y)
  abundance = sweep(Y, 1, nsize, "/")
  distance = vegan::vegdist(x = abundance, method = method)
  mds = ape::pcoa(distance)
  mds1 = mds$vectors[,1]
  mds2 = mds$vectors[,2]

  library(ggplot2)

  figure = list()
  #### With main_variable input
  if (!is.null(main_variable)){
    main_variable = as.factor(main_variable)
    point =data.frame(PC1 = mds1, PC2 = mds2, batch, main = main_variable)
    k = 1
    for (i in levels(point$main)){
      curve = data.frame()
      point_sub=point[point$main == i,]
      for (g in levels(point_sub$batch)){
        df = with(point_sub[point_sub$batch==g,], ellipse::ellipse(x = cor(PC1, PC2), scale = c(sd(PC1), sd(PC2)),
                                                           centre = c(mean(PC1), mean(PC2)),level = 0.8))
        df = data.frame(df,group = g)
        curve = rbind(curve,df)
      }
      g = ggplot(aes(x = PC1, y = PC2, color = batch),data = point_sub) +
          geom_point(size = 2) + theme(legend.position = "right") +
          ggtitle(paste("main_variable =", i)) +
          theme(plot.title = element_text(hjust = 0.5), title = element_text(size = 12)) +
          theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          theme(legend.text = element_text(size = 13)) +
          theme(legend.title = element_text(size = 13)) +
          theme(axis.title = element_text(size = 13)) +
          theme(axis.text = element_text(size = 12)) +
          geom_path(data = curve, aes(x = x, y = y, colour = group), size = 1, linetype = 1)
      figure[[k]] = g
      k = k + 1
    }
  }

  #### Without main_variable input
  if (is.null(main_variable)){
    point = data.frame(PC1 = mds1, PC2 = mds2, batch)
    curve = data.frame()
    for (g in levels(point$batch)){
      df = with(point[point$batch==g,], ellipse::ellipse(x = cor(PC1, PC2), scale = c(sd(PC1), sd(PC2)),
                                                centre = c(mean(PC1), mean(PC2)),level = 0.8))
      df = data.frame(df,group=g)
      curve = rbind(curve,df)
    }
    g = ggplot(aes(x = PC1, y = PC2, color = batch), data = point) +
        geom_point(size = 2) + theme(legend.position = "right") +
        theme(plot.title = element_text(hjust = 0.5),title = element_text(size = 12)) +
        theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(legend.text = element_text(size = 13)) +
        theme(legend.title = element_text(size = 13)) +
        theme(axis.title = element_text(size = 13)) +
        theme(axis.text = element_text(size = 12)) +
        geom_path(data = curve, aes(x = x, y = y, colour = group), size = 1, linetype = 1)
    figure = g
  }
  return(figure)
}




#' Threshold the posterior inclusion probability (PIP) through control Bayesian false discovery
#' rate (bFDR).
#' @param PIP_vec A vector contains the PIPs of parameters
#' @param alpha The level of the bFDR to need to control (default = 0.1)
#' @return The cutoff for PIPs to control the bFDR with the user defined value, alpha.
#' @examples
#' data(L_mean)
#' cutoff <- fdr_cut(L_mean, alpha = 0.1)
#' @export
fdr_cut = function(PIP_vec, alpha = 0.1){
  p = sort(1 - PIP_vec)
  thres = cumsum(p)/c(1:length(p))
  k = which.max(thres >= alpha)
  cutoff = 1 - p[k-1]
  return(cutoff)
}





#' Trace plot of BDMMA output
#' @param trace A data.frame named "trace" contained in the output of function BDMMA
#' @param param A character vector including the parameters' name for trace_plot
#' @param col A string defining the color of trace plot (default color is black)
#' @return The function returns a list containing plot objects of parameters' trace plot.
#' @examples
#' data(dat)
#' X <- dat$X
#' Y <- dat$Y
#' batch <- dat$batch
#' continuous <- dat$continuous
#' output <- BDMMA(X, Y, batch, continuous, burn_in = 3000, sample_period = 3000)
#' figure <- trace_plot(trace, param = c("alpha_1", "beta1_10"))
#' print(figure)
#' @export
#'
trace_plot = function(trace, param, col = "black"){
  for (i in param){
    if (!(i %in% names(trace))){
      stop(paste("The parameter", i, "is not in the list!"))
    }
  }
  param_trace = list()
  k = 1

  library(ggplot2)
  for (j in param){
    trace_j = trace[, names(trace) == j]
    trace_j = data.frame(iter = c(1: nrow(trace)), value = trace_j)
    g = ggplot(aes(x = iter, y = value), data = trace_j) +
      geom_line(colour = col) + ggtitle(paste("The trace plot of parameter", j)) +
      theme(plot.title = element_text(hjust = 0.5),title = element_text(size = 12)) +
      theme(legend.text = element_text(size = 13)) +
      theme(legend.title = element_text(size = 13)) +
      theme(axis.title = element_text(size = 13)) +
      theme(axis.text = element_text(size = 12))
    param_trace[[k]] = g
    k = k + 1
  }
  return(param_trace)
}










