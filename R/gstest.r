#' Generalized Two-Sample Kolmogorov-Smirnov Test for Right-Censored Data
#'
#' This function can perform the generalized two-sample Kolmogorov-Smirnov test
#' for right-censored data that proposed by Fleming et al. (1980).
#'
#' @name gstest
#' @rdname gstest
#' @param formula a formula object, which must have a \code{Surv} object as the
#' response on the left of the ~ operator and only one term on the right.
#' This test is a two-sample test, and this function only supports a two-level
#' covariate. This function treats 1 as an event and 0 as a censoring.
#' @param data a data frame in which to interpret the variables named in the
#' formula.
#' @return
#' \itemize{
#' \item \code{A}: Test statistic \eqn{A}; see Fleming et al. (1980).
#' \item \code{R}: Test statistic \eqn{R}; see Fleming et al. (1980).
#' \item \code{pval}: Two-sieded P-value.
#' }
#' @references
#' Fleming, D., O'Fallon, JR., O'Brien, P.C., and Harrington, DP. (1980).
#' Modified Kolmogorov-Smirnov test procedures with application to arbitrarily
#' right-censored data.
#' \emph{Biometrics}
#' \strong{36}(4): 607-625.
#' \url{https://doi.org/10.2307/2556114}
#' @examples
#' library(survival)
#'
#' # Example 1 of Fleming et al. (1980)
#' df1 <- data.frame(
#'   y = c(28,89,175,195,309,377,393,421,447,462,709,744,770,1106,1206,
#'         34,88,137,199,280,291,299,300,309,351,358,369,369,370,375,
#'         382,392,429,451,1119),
#'   event = c(1,1,1,1,1,0,0,0,0,1,0,0,0,0,0,
#'             1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,0,1,0),
#'   arm = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'           1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
#' )
#' # A = 1.834, R = 0.653, Two-sided P = 0.00223
#' gstest(Surv(y, event) ~ arm, df1)
#' @importFrom survival Surv
#' @importFrom stats pnorm model.frame model.response model.matrix
#' @export
gstest <- function(formula, data) {

  mf <- model.frame(formula, data)
  n <- nrow(mf)
  y <- model.response(mf)
  x <- model.matrix(terms(mf), mf)
  arm <- x[,2]
  armset <- unique(arm)

  if (dim(x)[2] > 2) {
    stop("Only a two-sample test is supported.")
  } else if (length(armset) != 2) {
    stop("Only a two-sample test is supported.")
  } else if (class(y) != "Surv") {
    stop("The outcome of the formula argument must be a \"Surv\" class object.")
  }

  time <- y[,1]
  event <- y[,2]
  ti <- unique(time[order(time)])

  # arm 1
  df1 <- data.frame(
    time = time[arm == armset[1]],
    event = event[arm == armset[1]]
  )
  n1 <- nrow(df1)
  df1 <- ties(df1)
  df1_stat <- stat_ab(df1, ti)
  names(df1_stat)[names(df1_stat) == "event"] <- "event1"
  names(df1_stat)[names(df1_stat) == "censor"] <- "censor1"
  names(df1_stat)[names(df1_stat) == "n"] <- "n1"
  names(df1_stat)[names(df1_stat) == "d"] <- "d1"
  names(df1_stat)[names(df1_stat) == "beta"] <- "beta1"
  names(df1_stat)[names(df1_stat) == "alpha"] <- "alpha1"

  # arm 2
  df2 <- data.frame(
    time = time[arm == armset[2]],
    event = event[arm == armset[2]]
  )
  n2 <- nrow(df2)
  df2 <- ties(df2)
  df2_stat <- stat_ab(df2, ti)
  names(df2_stat)[names(df2_stat) == "event"] <- "event2"
  names(df2_stat)[names(df2_stat) == "censor"] <- "censor2"
  names(df2_stat)[names(df2_stat) == "n"] <- "n2"
  names(df2_stat)[names(df2_stat) == "d"] <- "d2"
  names(df2_stat)[names(df2_stat) == "beta"] <- "beta2"
  names(df2_stat)[names(df2_stat) == "alpha"] <- "alpha2"

  df3 <- merge(df1_stat, df2_stat, by = "time")

  df3$eta <- 1/sqrt(1.0/(n1*exp(-df3$alpha1)) + 1.0/(n2*exp(-df3$alpha2)))

  df3$u <- 0
  u <- 0
  for(i in 1:nrow(df3)) {
    if (df3$event1[i] >= 1 & df3$event2[i] >= 1) {
      l1 <- 0
      l2 <- 0
      for (j in 1:df3$event1[i]) {
        l1 <- l1 + 1.0/(df3$n1[i] - (j - 1))
      }
      for (j in 1:df3$event2[i]) {
        l2 <- l2 + 1.0/(df3$n2[i] - (j - 1))
      }
      u <- u + df3$eta[i]*(l1 - l2)
    } else if (df3$event1[i] >= 1) {
      l1 <- 0
      for (j in 1:df3$event1[i]) {
        l1 <- l1 + 1.0/(df3$n1[i] - (j - 1))
      }
      u <- u + df3$eta[i]*(l1 - 0)
    } else if (df3$event2[i] >= 1) {
      l2 <- 0
      for (j in 1:df3$event2[i]) {
        l2 <- l2 + 1.0/(df3$n2[i] - (j - 1))
      }
      u <- u + df3$eta[i]*(0 - l2)
    }
    df3$u[i] <- u
  }
  df3$ynt <- (exp(-df3$beta1) + exp(-df3$beta2))*df3$u*0.5

  J <- max((1:nrow(df3))[df3$event1 >= 1], (1:nrow(df3))[df3$event2 >= 1])
  R <- 1 - (exp(-df3$beta1[J]) + exp(-df3$beta2[J]))*0.5

  # one-sided test
  # V <- max(0, df3$ynt)
  # pvr <- 1 - pnorm(V/sqrt(R - R^2)) + pnorm(V*(2*R - 1)/sqrt(R - R^2))*exp(-2*V^2)

  # two-sided
  A <- max(abs(df3$ynt))
  par <- 2*(1 - pnorm(A/sqrt(R - R^2)) + pnorm(A*(2*R - 1)/sqrt(R - R^2))*exp(-2*A^2))

  data.frame(
    A = A,
    R = R,
    pval = par
  )

}

ties <- function(df) {
  tbl <- table(df$time, df$event)
  tbl <- as.data.frame(tbl)
  tbl0 <- tbl[tbl$Var2 == 0, c("Var1", "Freq")]
  tbl1 <- tbl[tbl$Var2 == 1, c("Var1", "Freq")]
  names(tbl0)[names(tbl0) == "Var1"] <- "time"
  names(tbl0)[names(tbl0) == "Freq"] <- "censor"
  names(tbl1)[names(tbl1) == "Var1"] <- "time"
  names(tbl1)[names(tbl1) == "Freq"] <- "event"
  tbl2 <- merge(tbl1, tbl0, by = "time")
  tbl2 <- tbl2[order(tbl2$time),]
  tbl2
}

stat_ab <- function(df0, ti) {
  n <- sum(df0$event) + sum(df0$censor)
  b <- 0
  a <- 0
  df1 <- data.frame(time = ti, n = NA, d = NA, beta = NA, alpha = NA)
  df <- merge(df1, df0, by = "time", all = TRUE)
  df$event[is.na(df$event)] <- -1
  df$censor[is.na(df$censor)] <- -1
  df[is.na(df)] <- 0
  for (i in 1:nrow(df)) {
    if (df$event[i] == -1) {
      df$n[i] <- n
      df$d[i] <- 0
      df$beta[i] <- b
      df$alpha[i] <- a
    } else {
      if (df$event[i] > 0) {
        for (j in 1:df$event[i]) {
          b <- b + 1.0/(n - (j - 1))
        }
      } else if (df$censor[i] > 0) {
        for (j in 1:df$censor[i]) {
          a <- a + 1.0/(n - (j - 1))
        }
      }
      df$n[i] <- n
      df$d[i] <- df$event[i]
      df$beta[i] <- b
      df$alpha[i] <- a
      n <- n - df$event[i] - df$censor[i]
    }
  }
  df
}
