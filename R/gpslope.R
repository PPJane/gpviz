#' A Gaussian process prediction reference on slopes
#'
#' @export
#' @param data A list of input values
#' @param x_ name of x variable, should match with the data frame variable name
#' @param y_ name of y variable, should match with the data frame variable name
#' @param trainingpect the percent of training instances, the default value is 10\% of the whole sample, avoid high number of training instances when you want the Stan program to be done faster
#' @param seed a random generator seed to designate such that your results will be the same across different runs
#' @param dn The number you want to sample for density calculation for the first derivative's distribution for each x*, the default value is 64
#' @param iter Number of iterations you would like Stan to run
#' @param chains Number of chains you would like Stan to utilize
#' @param conf Confidence level resulting in the interval displayed for the first derivative of the gaussian process prediction
#' @param condition1 The first condition you want the first derivative to satisfy, such as ">0", the default value is NULL
#' @param condition2 The second condition you want the first derivative to satisfy, such as "<2", the default value is NULL
#' @param rho_alpha Prior distribution of rho is set to inverse gamma distribution with rho ~ inv_gamma(rho_alpha, rho_beta), default value is inv_gamma(4, 4)
#' @param rho_beta Prior distribution of rho is set to inverse gamma distribution with rho ~ inv_gamma(rho_alpha, rho_beta), default value is inv_gamma(4, 4)
#' @param alpha_mean Prior distribution of alpha is set to normal distribution with alpha ~ normal(alpha_mean, alpha_sd), default value is normal(0, 1)
#' @param alpha_sd Prior distribution of alpha is set to normal distribution with alpha ~ normal(alpha_mean, alpha_sd), default value is normal(0, 1)
#' @param sigma_mean Prior distribution of sigma is set to normal distribution with sigma ~ normal(sigma_mean, sigma_sd), default value is normal(0, 1)
#' @param sigma_sd Prior distribution of sigma is set to normal distribution with sigma ~ normal(sigma_mean, sigma_sd), default value is normal(0, 1)
#' @return A list consiting of the Stan fit object, three ggplot plot objects and the user-specified argument values
#' @examples resultlst <- gpviz(data = hsbReport_cleaned, x_ = 'deposit_asset_ratio', y_ = 'loss_asset_ratio', seed = 123456, iter = 2000, conf = 0.95, condition1 = '>0', condition2 = NULL)

gpslope <- function(data, x_, y_, trainingpect = 0.1, seed = 123,
                   iter = 100, chains = 4,
                   conf = 0.5, condition1 = NULL, condition2 = NULL,
                   dn = 64,
                   rho_alpha = 4, rho_beta = 4, alpha_mean = 0, alpha_sd = 1, sigma_mean = 0, sigma_sd = 1){
  arg <<- as.list(match.call())

  data <- data %>% select_at(c(x_, y_)) %>% arrange_at(x_)

  N1 <- floor(trainingpect * nrow(data))

  set.seed(seed)
  train_idx <- sort(sample(1:nrow(data), size = N1))
  train_data <- data[train_idx, ]
  test_data <- data

  x1 <- train_data %>% select(x_) %>% unlist()
  y1 <- train_data %>% select(y_) %>% unlist()

  x2o <- data %>% select(x_) %>% unlist()
  y2o <- data %>% select(y_) %>% unlist()

  x2 <- union(x2o, seq(from = range(x2o)[1], to = range(x2o)[2], length.out = 100))
  x2 <- sort(x2)
  N2 <- length(x2)

  pred_data <- list(N1 = N1, N2 = N2, x1 = x1, y1 = y1, x2 = x2, y2o = y2o, x2o = x2o, x_ = x_, y_ = y_, rho_alpha, rho_beta, alpha_mean, alpha_sd, sigma_mean, sigma_sd)

  pred_fit <- rstan::stan(file = '../src/stan_files/predict_gauss.stan', data = pred_data, iter = iter, chains = chains)
  print(pred_fit, pars = c('alpha', 'rho', 'sigma'))

  f2 <- rstan::extract(pred_fit)$f2
  f2_prime <- rstan::extract(pred_fit)$f2_prime

  confInt <- function(x){
    quantile(x, probs = c((1 - conf)/2, (1 + conf)/2))
  }
  confIntDf <- apply(f2_prime, 2, confInt)
  confIntDf <- data.frame(cbind(x = x2, t(confIntDf)))
  colnames(confIntDf) <- c('x', 'low', 'high')

  if (is.null(condition1)){
    print("Warning: You did not specify any requirement for the first derivative, such as '> 0'.")
    roi <- data.frame()
  } else if (is.null(condition2)){
    roi <- confIntDf %>% filter(eval(parse(text = paste('low', condition1, '& high', condition1))))
    roix <- roi %>% select(x)
    condition <- TeX(paste("$\\frac{ \\partial f(x^{*})}{ \\partial x^{*}}$", condition1))
  } else {
    roi <- confIntDf %>% filter(eval(parse(text = paste('low', condition1, '& high', condition1, '& low', condition2, '& high', condition2))))
    roix <- roi %>% select(x)
    condition <- TeX(paste("$\\frac{ \\partial f(x^{*})}{ \\partial x^{*}}$", condition1, " & $\\frac{ \\partial f(x^{*})}{ \\partial x^{*}}$", condition2))
  }

  p1 <- ggplot2::ggplot() +
    geom_point(data = data.frame(x = x2o, y = y2o), aes(x, y, color = 'test'), alpha = 0.1) +
    geom_point(data = data.frame(x = x1, y = y1), aes(x, y, color = 'training'), alpha = 0.8) +
    geom_point(data = data.frame(x = x2, f2 = colMeans(f2)), aes(x2, f2, color = 'prediction'), alpha = 0.9) +
    scale_color_manual(name = '',
                       labels = list('test' = 'Test \nData',
                                     'training' = 'Training \nData',
                                     'prediction' = 'Prediction \nExpectation',
                                     'realization' = 'Randomly \nPicked \nRealization'),
                       values = c('test' = 'blue',
                                  'training' = 'purple',
                                  'prediction' = 'pink',
                                  'realization' = 'grey')) +
    ylab(gsub('_', '_\n', y_)) +
    ggtitle(paste0("Gaussian Process and the Confidence Interval of \nIts Prediction's First Derivative at Confidence Level of ",
                  sprintf('%i%%', conf * 100),
                  '\nFor Data Set "',
                  # arg$data)) +
                  ifelse(exists("argVal"), argVal$data, arg$data), '"')) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16, angle = 0, vjust = 0.5),
          title = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 0.5),
          legend.position = 'right',
          legend.text = element_text(size = 12),
          legend.key = element_rect(size = 5),
          legend.key.size = unit(3, 'lines'))

  p2 <- ggplot2::ggplot(data = confIntDf) +
    geom_ribbon(aes(x, ymin = low, ymax = high), fill = 'grey', alpha = 0.4) +
    scale_color_manual(name = '',
                       labels = list('low' = 'Lower \nConfidence \nInterval',
                                     'high' = 'Upper \nConfidence \nInterval',
                                     'roi' = 'Region of \nInterest'),
                       values = c('low' = 'yellow',
                                  'high' = 'orange',
                                  'roi' = 'pink')) +
    xlab(x_) +
    ylab(TeX("$\\frac{ \\partial f(x^{*})}{ \\partial x^{*}}$")) +
    theme(axis.title.x = element_text(size = 16, vjust = 0.5),
          axis.title.y = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 0.5),
          legend.position = 'right',
          legend.text = element_text(size = 12),
          legend.key = element_rect(size = 5),
          legend.key.size = unit(3, 'lines'))

  if (nrow(roi) > 0){
    p2 <- p2 +
      geom_segment(data = roi, aes(x = x, xend = x, y = low, yend = high, color = 'roi'), size = 2) +
      geom_point(data = data.frame(x = roix,
                                   y = ggplot_build(p2)$layout$panel_ranges[[1]]$y.range[1]),
                 aes(x, y, color = 'roi')) +
      labs(subtitle = condition)
  }

  p2 <- p2 +
    stat_smooth(aes(x, y = low, color = "low"), method = 'lm', formula = y ~ poly(x,15), size = 1.5, se = FALSE, alpha = 0.01) +
    stat_smooth(aes(x, y = high, color = 'high'), method = 'lm', formula = y ~ poly(x,15), size = 1.5, se = FALSE, alpha = 0.01)

  for (i in sample(1:nrow(f2), 10)) {
    p1 <- p1 +
      geom_point(data = data.frame(x = x2, y = unlist(f2[i,])), aes(x, y, color = 'realization'), alpha = 0.1)
  }

  denDf <- data.frame(NA, nrow = N2 * dn, ncol = 3)
  for (i in 1:N2){
    tmp <- density(f2_prime[,i], n = dn)
    denDf[((i-1)*dn+1):(i*dn),] <- data.frame(x = x2[i], y = tmp$x, den = tmp$y)
  }
  colnames(denDf) <- c('x', 'f2p', 'den')

  palette.breaks <- quantile(denDf$den, probs = seq(0, 1, by = 0.1))
  palette.limits <- quantile(denDf$den, probs = c(0, 0.975))
  palette.color  <- colorRampPalette(c("#DADADA", "#FFFFBF","#FC8D59"))(length(palette.breaks) - 1)

  p3 <- ggplot2::ggplot(data = denDf) +
    geom_point(aes(x, f2p, colour = den)) +
    scale_color_gradientn(colours = palette.color,
                          breaks = palette.breaks,
                          limits = palette.limits,
                          na.value = 'red',
                          guide = guide_colorbar(title =
                                                   TeX('$P \\left( \\frac{ \\partial f(x^{*})}{ \\partial x^{*}} \\right)$ in Percentile'),
                                                 barheight = 10,
                                                 barwidth = 1,
                                                 nbin = 10)) +
    xlab(x_) +
    ylab(TeX("$\\frac{ \\partial f(x^{*})}{ \\partial x^{*}}$")) +
    theme(axis.title.x = element_text(size = 16, vjust = 0.5),
          axis.title.y = element_text(size = 16, angle = 0, vjust = 0.5),
          legend.position = 'right',
          legend.text = element_text(size = 12),
          legend.key = element_rect(size = 5),
          legend.key.size = unit(3, 'lines'))

  return(list(pred_fit, p1, p2, p3, arg))
}
