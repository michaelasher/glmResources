#' Determine the minimum detectable effect size
#'
#' Given desired power, degrees of freedom, and alpha, this function returns the
#' minimum detectable effect size as an f-squared statistic.
#'
#' @param desired_power Defaults to .8, but this can be changed to test for the
#' minimum detectable effect size for other levels of power.
#' @param num_df Numerator degrees of freedom for the test.
#' @param den_df Denominator degrees of freedom for the test.
#' @param alpha Alpha level for the test.
#' @return Returns an f-squared statistic representing the minimum detectable effect
#' size with the desired power, degrees of freedom, and alpha.
mdes <- function(desired_power = .8, num_df, den_df, alpha = .05) {

  # overshoot, nearest tenth
  power <- 0
  fsq <- 0
  while (power < desired_power) {
    fsq <- fsq + .01
    power <- pwr::pwr.f2.test(u = num_df, v = den_df, f2 = fsq, sig.level = alpha)$power
  }

  # overshoot, nearest thousandth
  fsq <- fsq - .01
  power <- 0
  while (power < desired_power) {
    fsq <- fsq + .0001
    power <- pwr::pwr.f2.test(u = num_df, v = den_df, f2 = fsq, sig.level = alpha)$power
  }

  # overshoot, nearest ten thousandth
  fsq <- fsq - .0001
  power <- 0
  while (power < desired_power) {
    fsq <- fsq + .00001
    power <- pwr::pwr.f2.test(u = num_df, v = den_df, f2 = fsq, sig.level = alpha)$power
  }

  return(fsq)

}


#' Calculate power and for GLM tests
#'
#' Wrapper to calculate power for tests of parameter estimates or full model in
#' GLM based on Cohen's tables and using pwr.f2.test in pwr package. Inspired by
#' modelPower in John Curtin's lmSupport package. Allows the use
#' of partial eta squared or delta R2 rather than just f2 as effect size. If
#' you provide power, it returns N, if you provide N, it returns power. If you
#' provide N and power, it returns the minimum detectable effect size given the
#' specified N and power. You must alwasy specify an effect size as either
#' f2, partial eta2, or delta R2 with model R2. You must also specify the
#' number of parameters in the compact (pc) and  augmented (pa) for the model
#' comparison that will test the effect.
#'
#' @param pc Number of parameters in the compact model; i.e., intercept + all
#' parameters excluding the effect of interest; This is the numerator df of the
#' F test for the effect.
#' @param pa Number of parameters in the augmented model; i.e., the intercept
#' and all parameters including the effect of interest.
#' @param n Sample size.
#' @param alpha Alpha for statistical test.
#' @param power Power for statistical test.
#' @param f2 Effect size.
#' @param peta2 = Partial eta2 effect size.
#' @param dr2 Delta r2 effect; if provided must also specify r2.
#' @param r2 Model r2, only needed if using delta r2.
#' @return Returns a list with n, power, possibly minimum detectable effect size.
#' @examples
#' # return the minimum detectable effect size with 200 participants, power of
#' # .8 (the default), pa of 5, and pc of 4:
#' power_analysis(pa = 5, pc = 4, n = 200)
#'
#' # return the number of participants needed for 70% power given
#' # pa of 3, pc of 2, and peta2 of .01
#' power_analysis(pa = 3, pc = 2, peta2 = .01, power = .7)
#'
#' # return the power of a study with peta2 of .02 and 50 participants, with
#' # pa of 5 and pc of 4.
#' power_analysis(pa = 5, pc = 4, peta2 = .02, n = 50)
power_analysis <- function(pc = NULL, pa = NULL, n = NULL, alpha = 0.05, power = NULL,
                           f2 = NULL, peta2 = NULL, dr2 = NULL, r2 = NULL) {
  if (is.null(pa) | is.null(pc)) {
    stop("Must provide pa and pc")
  }

  u <- pa - pc
  mdes_peta2 <- NULL
  mdes_f2 <- NULL

  nEffs <- 0
  if (!is.null(f2)) {
    nEffs <- nEffs + 1
  }
  if (!is.null(peta2)) {
    f2 <- peta2 / (1 - peta2)
    nEffs <- nEffs + 1
  }
  if (!is.null(dr2) & !is.null(r2)) {
    f2 <- dr2 / (1 - r2)
    nEffs <- nEffs + 1
  }
  if (nEffs > 1) {
    stop("Must not specify more than one of the following: f2, peta2, or both dr2 and r2")
  }

  if (!is.null(n)) {
    v <- n - pa
  }
  else {
    v <- NULL
  }

  if (!is.null(pa) & !is.null(pc) & !is.null(n)) {
    mdes_power <- power
    if(is.null(mdes_power)) mdes_power <- .8
    mdes_f2 <- mdes(desired_power = mdes_power, u, v, alpha)
    mdes_peta2 <- mdes_f2 / (1 + mdes_f2)
  }

  if (nEffs != 0) {
    results <- pwr::pwr.f2.test(
      u = u, v = v, f2 = f2, sig.level = alpha,
      power = power
    )

    output <- list(
      n = pa + results$v,
      peta2 = results$f2 / (1 + results$f2),
      f2 = results$f2,
      power = results$power,
      mdes_peta2 = mdes_peta2,
      mdes_f2 = mdes_f2
    )
  }

  if (nEffs == 0) {
    print("2")
    output <- list(
      n = n,
      peta2 = NULL,
      f2 = NULL,
      power = power,
      mdes_peta2 = mdes_peta2,
      mdes_f2 = mdes_f2
    )
  }

  return(output)
}
