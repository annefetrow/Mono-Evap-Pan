#' York regression
#' 
#' Adapted from the IsoplotR::york and geostats::york implementations which are direct translations of York, D., Evensen, N.M., Martı́nez, M.L., De Basabe Delgado, J. (2004) Unified equations for the slope, intercept, and standard errors of the best straight line. American Journal of Physics 72, 367–375.
#' Additional useful literature:
#'  - Wehr, R., Saleska, S.R. (2017) The long-solved problem of the best-fit straight line: Application to isotopic mixing lines. Biogeosciences 14, 17–29.
#'  - York, D. (1968) Least squares fitting of a straight line with correlated errors. Earth and Planetary Science Letters 5, 320–324.
#'  - see https://github.com/JENScoding/yorkregression for other potential output and return value that actually fits an R type model (and might be used with broom functions)
#' 
#' Features:
#' - input as vectors instead of alpha matrix and output as tibble instead of list to be more tidyverse-friendly, i.e. easily used with group_by + summarize + unnest of the result
#' 
#' @param X/Y x/y coordinates of the observed data points
#' @param sigma_X/Y errors in the x/y coordinates
#' @param omega_X/Y weights of each data point, default is 1/sigma^2 (see step #2 in York et al. 2004)
#' @param r_X/Y the correlation coefficient between the errors in X and Y for each data point
#' @param maxiter alpha positive integer specifying the maximum number of iterations allowed.
#' @param tol alpha positive numeric value specifying the desired tolerance for successive estimates of b (see step #6 in York et al. 2004)
#' @return 
york_regression <- function (X, Y, sigma_X, sigma_Y, omega_X = 1/sigma_X^2, omega_Y = 1/sigma_Y^2, r_XY = 0, maxiter = 50, tol = 1e-15) {
  
  # make alpha simple data frame to also get dimensionality checks
  ds <- dplyr::tibble(
    X, Y, omega_X, omega_Y, r_XY
  ) |> dplyr::filter(!is.na(X), !is.na(Y), !is.na(omega_X), !is.na(omega_Y))
  
  # safety check if enough data
  if (nrow(ds) < 2) {
    warning("not enough data for fit", immediate. = TRUE, call. = FALSE)
    return(tibble::tibble(chisq_df = nrow(ds) - 2L))
  }
  
  # initial guess for alpha (intercept) and b (slope) from simple linear model
  # see step #1 in York et al. 2004
  ab <- stats::lm(Y ~ X, data = ds)$coefficients
  b <- ab[2]
  
  # safety check
  if (any(is.na(ab)) || b == 0) {
    warning("no straight line fit possible", call. = FALSE)
    return(tibble::tibble(a = ab[1], b = b, chisq_df = nrow(ds) - 2L))
  }
  
  # calculate constat fitting parameters
  ds <- ds |> dplyr::mutate(
    omega_XY = omega_X * omega_Y,
    alpha = sqrt(omega_XY)
  )
  
  # iteratively find best b (steps #3-6 in York et al. 2004)
  for (i in 1:as.integer(maxiter)) {
    # calculate W (step #3 in York et al. 2004)
    W <- ds$omega_XY/(ds$omega_X + b^2 * ds$omega_Y - 2 * b * ds$r_XY * ds$alpha)
    # calculate derived quantities (step #4 in York et al. 2004)
    Xbar <- sum(W * ds$X, na.rm = TRUE)/sum(W, na.rm = TRUE)
    Ybar <- sum(W * ds$Y, na.rm = TRUE)/sum(W, na.rm = TRUE)
    U <- ds$X - Xbar
    V <- ds$Y - Ybar
    beta <- W * (U/ds$omega_Y + b * V/ds$omega_X - (b * U + V) * ds$r_XY/ds$alpha)
    # calculate new estimate of b (step #5 in York et al. 2004, eq. 13b)
    new_b <- sum(W * beta * V, na.rm = TRUE)/sum(W * beta * U, na.rm = TRUE)
    # check convergence (step #6 in York et al. 2004)
    converged <- abs(new_b/b - 1) < tol
    b <- new_b
    if (converged) break
  }
  
  # calculate final value of a (step #7 in York et al. 2004, eq 13a)
  a <- Ybar - b * Xbar
  
  # calculate adjusted x/y values (step #8 in York et al. 2004)
  xadj <- Xbar + beta
  
  # calculate final xbar and u (step #9 in York et al. 2004)
  xbar <- sum(W * xadj, na.rm = TRUE)/sum(W, na.rm = TRUE)
  u <- xadj - xbar
  
  # calculate error estimates of a and b (step #10 in York et al. 2004, et. 13c&d)
  sigma_b <- sqrt(1/sum(W * u^2, na.rm = TRUE))
  sigma_a <- sqrt(1/sum(W, na.rm = TRUE) + (xbar * sigma_b)^2)
  
  # goodness of fit (see York et al. 2004, section V)
  S <- sum(W * (ds$Y - b * ds$X - a)^2) # this is a chi^2 statistic
  chisq_df <- nrow(ds) - 2L # degrees of freedom

  # other stats
  # reduced chi^2 statistic (mean square weighted deviation in geochron community)
  # see https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
  # note that this is reported as mswd in IsoplotR::york()
  red_chisq <- S/chisq_df
  # the standard error for the chi^2 statistic
  sigma_red_chisq <- sqrt(2/chisq_df)
  # p_value 
  # note that this is the p-value for the H0 hypothesis that a linear regression is NOT a good fit
  # if this is high you CANNOT reject the H0 hypothesis (which means the linear reg is indeed NOT a good fit)
  # for the opposite (H0: linear regression IS a good fit), calculate 1 - p_chisq_p_val
  # the latter (1-p_val) is how the p-value is report in IsoplotR::york()!
  york_pval <- stats::pchisq(S, chisq_df) 
  
  # assemble return values
  retval <- tibble::tibble(
    # range
    xmin = min(ds$X), xmax = max(ds$X),
    ymin = min(ds$Y), ymax = max(ds$Y),
    # parameters
    a, b, sigma_a, sigma_b,
    # covariance
    cov_ab = -xbar * sigma_b^2,
    # correlation coefficient
    r_ab = -xbar * sigma_b/sigma_a,
    # chisq stats
    chisq = S, chisq_df, red_chisq, sigma_red_chisq, york_pval
  )
  
  return(retval)
}

#' generate a york regression line based on 
#' @param ds data frame with york regression parameter columns 'a', 'b', 'cov_ab', 'sigma_a', 'sigma_b', 'xmin', 'xmax'
#' @param expand how far to expand the line beyond the data's xmin/xmax (see [dplyr::expansion]), default is 10%
#' @param xmin/xmax define a specific xmin/xmax to model, if set @expand has no effect
generate_york_regression_line <- function(ds, expand = ggplot2::expansion(mult = 0.1), xmin = NULL, xmax = NULL, n = 10L) {
  
  expand_range <- function(min, max, mult, add, sign) {
    if (sign < 0) min - (max - min) * mult - add
    else max + (max - min) * mult + add
  }
  calc_x_from_y <- function(y, a, b) (y - a)/b
  
  # not used yet
  # calc_x_with_ci <- function(y, a, b) {
  #   (cov_ab + b * y - a * b + sqrt( (a * b - b * y - cov_ab)^2 - (b^2 - sb^2) * (a^2 - 2 * a * y - sa^2 + y2) ) ) / (b^2 - sb^2)
  # }
  
  ds |>
    dplyr::mutate(
      # xmin
      yrl_xmin_x = expand_range(xmin, xmax, expand[1], expand[2], -1),
      yrl_xmin_y = expand_range(calc_x_from_y(ymin, a, b), calc_x_from_y(ymax, a, b), expand[1], expand[2], -1),
      # xmax
      yrl_xmax_x = expand_range(xmin, xmax, expand[3], expand[4], +1),
      yrl_xmax_y = expand_range(calc_x_from_y(ymin, a, b), calc_x_from_y(ymax, a, b), expand[3], expand[4], +1),
      # is this always the right choice?
      yrl_xmin = if(!is.null(!!xmin)) !!xmin else yrl_xmin_y, # is this always the right choice?
      yrl_xmax = if(!is.null(!!xmax)) !!xmax else yrl_xmax_y
      #yrl_xmin = ifelse(yrl_xmin1 > yrl_xmin2, yrl_xmin1, yrl_xmin2),
      #yrl_xmax = ifelse(yrl_xmax_x < yrl_xmax_y, yrl_xmax_x, yrl_xmax_y)
    ) |>
    dplyr::select(-starts_with("yrl_xmin_"), -starts_with("yrl_xmax_")) |>
    dplyr::cross_join(data.frame(yrl_i = 1:n)) |>
    dplyr::mutate(
      # y = a + x * b --> x = (y - a) / b
      # sigma_y = sqrt(sa^2 + sb^2 * x^2 + 2 * x * cov_ab)
      # sigma_x = (y - a) / b * sqrt(sa^2/(y-a)^2 + sb^2 / b^2 + 2 * 1/b * 1/(y-a) * cov)
      # ymax = a + x * b + sqrt(sa^2 + sb^2 * x^2 + 2 * x * cov)
      # ymin = a + x * b - sqrt(sa^2 + sb^2 * x^2 + 2 * x * cov) 
      # --> x = (cov + b * y - a * b +/- sqrt( (a*b - b*y - cov)^2 - (b^2 - sb^2) * (a^2 - 2 * a * y - sa^2 + y^2) ) ) / (b^2 - sb^2)
      yrl_x = yrl_xmin + (yrl_i - 1L) * (yrl_xmax - yrl_xmin) / (!!n - 1L),
      yrl_y = a + b * yrl_x,
      yrl_sigma_y = sqrt(sigma_a^2 + (sigma_b * yrl_x)^2 + 2 * yrl_x * cov_ab),
      yrl_sigma_x = abs((yrl_y - a) / b) * sqrt( (sigma_a/(yrl_y - a))^2 + (sigma_b/b)^2 + 2 * 1/b * 1/(yrl_y - a) * cov_ab)
    )
}
