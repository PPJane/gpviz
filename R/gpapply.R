#' A Gaussian process prediction reference function applied over a list of features
#'
#' @export
#' @param data A list of input values
#' @param x_lst a list of names of x variables to be examed, they should match with the data frame variable names
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
#' @return Pdf plots describing Gaussian process prediction and relationship intensity in uncertainty interval and density estimation
#' @examples gpapply(data = winequality_red, x_lst = wine_red_x_lst, y_ = 'quality', iter = 1000, condition1 = '>0')

gpapply <- function(data, x_lst, y_, trainingpect = 0.1, seed = 123,
                    iter = 100, chains = 4,
                    conf = 0.5, condition1 = NULL, condition2 = NULL,
                    dn = 64,
                    rho_alpha = 4, rho_beta = 4, alpha_mean = 0, alpha_sd = 1, sigma_mean = 0, sigma_sd = 1){
  argVal <<- as.list(sys.call())

  foo <- function(feature){
    tmplist <- gpviz(data, feature, y_, trainingpect, seed,
                     iter, chains,
                     conf, condition1, condition2,
                     dn,
                     rho_alpha, rho_beta, alpha_mean, alpha_sd, sigma_mean, sigma_sd)

    gpPlot <- tmplist[[2]]
    gpSlopeCI <- tmplist[[3]]
    gpSlopeDen <- tmplist[[4]]

    pdf(file = paste0("./gp_CI_", argVal$data, "_",
                      ifelse(is.numeric(argVal$conf), argVal$conf, formals(gpviz)$conf), "_",
                      feature, "_", argVal$y_, ".pdf"),
        height = 20, width = 18, paper = "letter")
    grid::grid.newpage()
    grid::grid.draw(rbind(ggplotGrob(gpPlot), ggplotGrob(gpSlopeCI), size = "last"))
    dev.off()


    pdf(file = paste0("./gp_density_", argVal$data, "_", feature, "_", argVal$y_, ".pdf"),
        height = 20, width = 18, paper = "letter")
    grid::grid.newpage()
    grid::grid.draw(rbind(ggplotGrob(gpPlot), ggplotGrob(gpSlopeDen), size = "last"))
    dev.off()
  }

  lapply(x_lst, foo)
}
