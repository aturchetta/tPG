
#' Average recruitment curve
#'
#' @description This function generates the average recruitment curve according to either the pdf or cdf of a
#' Gamma random variable to be used in the \code{\link{sim_rec}} function to simulate recruitment data.
#' @param p1 Shape parameter.
#' @param p2 Rate parameter.
#' @param t Grid of time points.
#' @param c2 Multiplicative constant.
#' @param c1 Shifting constant.
#' @param type Either the pdf of a Gamma random variable (default) or the cdf.
#'
#' @importFrom stats pgamma rgamma rpois qnorm
#' @importFrom graphics lines points
#' @import splines2 nloptr
#'
#' @return A vector with the average recruitment intensity at each time point.
#' @export
#'
#' @examples
lambda <- function(p1, p2, t, c2 ,c1, type){

  if(type=="cdf"){

    c2*pgamma(t, shape = p1 , rate = p2) + c1

  }else{

    c2*(p2^p1)*t^{p1-1}*exp(-p2*t)/gamma(p1) + c1

  }


}






#' Data generating function
#'
#' @param alpha Shape parameter of the Gamma distribution.
#' @param p1 Shape parameter of the average recruitment intensity.
#' @param p2 Rate parameter of the average recruitment intensity.
#' @param c2 Multiplicative constant of the average recruitment intensity.
#' @param c1 Shifting constant of the average recruitment intensity.
#' @param type Form of the average recruitment intensity. It can be either the pdf of a Gamma random variable (default) or the cdf.
#' @param t_p Plateau point.
#' @param Tmax Maximum value for the grid of time points.
#' @param s_time Vector of the centers' initiation times. Must be of length \code{C_act} if \code{random_init = F}.
#' @param C_act Number of active centers.
#' @param random_init Logical; if TRUE (default), the centers' initiation times are selected randomly from the values in \code{s_time}.
#'
#' @return A list with two elements: the vector of initiation times and a list where each element contains the recruitment data of one center
#' formatted as a data frame with two columns; the first column is the grid of recruitment days for the concerned center,
#' and the second column indicates the respective number of enrollments on that day.
#' @export
#'
#' @examples \dontrun{sim_data <- sim_rec(alpha = 1, p1 = 10, p2 = 0.15, c2 = 0.3, c1 = 0.1,
#' type = "cdf", t_p = 150, Tmax = 200, s_time = 1:100, C_act = 30, random_init = TRUE)}
sim_rec <- function(alpha, p1, p2,  c2 , c1 , type , t_p , Tmax , s_time , C_act, random_init = TRUE){

  start <- c(1,sort(sample(s_time, size = C_act-1, replace = random_init)))

  t <- 1:Tmax
  t.seq <- 1:t_p

  cent <- list()

  centpred <- matrix(NA, nrow = length(t), ncol = C_act)

  for (j in 1:C_act) {

    lam_cen <- rgamma(length(t.seq),alpha,alpha/lambda(p1 = p1, p2 = p2, t = t.seq,
                                                       c2 = c2, c1 = c1, type = type))
    lam_cen <- c(lam_cen , rep(lam_cen[length(t.seq)],  length(t)-length(t.seq) ))
    rec_cen <- rpois(length(t), lambda = lam_cen)

    cent[[j]] <- data.frame(t = t[1:(Tmax-start[j]+1)],
                            count = rec_cen[1:(Tmax-start[j]+1)])
  }

  return(list(data = cent , start = start))
}







#' Estimating function of the tPG model parameters
#'
#' @param data Recruitment data of the study formatted as follows: a list where each element contains the recruitment data of one center formatted as a data frame with two columns;
#' the first column is the grid of recruitment days (or another time unit) for the concerned center, and the second column indicates the respective number of enrollments on that day.
#' @param t_int Interim time of the analysis with respect to the start of the first center.
#' @param tp_start Starting values for the plateau point to be used in the optimization function. The likelihood is optimized for all these starting values and the the estimated parameters
#' corresponding to the highest likelihood are used.
#' @param degree Candidate degree values for the polynomial.
#' @param den Candidate internal knot placement, expressed as a function of the estimated plateau point, i.e., t_p/\code{den}
#' @param max_eval Maximum number of iterations of the optimization function. Default is 100000.
#' @param tol_rel Tolerance of the optimization function. Default is 1.e-8.
#' @param print_models Logical. Print iterations across internal knots and polynomial degrees candidates.
#' @param low_tp Lower bound for the plateau point to be used in the optimization function.
#'
#' @return This function returns multiple objects which include the estimates of the parameters under each model and the best-fitting model selected via BIC.
#' @export
#'
#' @examples \dontrun{sim_data <- sim_rec(alpha = 1, p1 = 10, p2 = 0.15, c2 = 0.3, c1 = 0.1,
#' type = "cdf", t_p = 150, Tmax = 200, s_time = 1:100, C_act = 30, random_init = TRUE)
#'
#' est <- tPG_est(data = sim_data$data, t_int = 200 , tp_start = c(20, 60, 100, 140), degree = c(2, 3), den = c(NA, 2, 3) )}
tPG_est <- function(data, t_int, tp_start, degree, den, max_eval=100000, tol_rel=1.e-8, print_models = TRUE , low_tp = 6){




  t <- 1:t_int


  results <- rep(list(list()),length(degree))
  bic_mat <- matrix(NA, nrow = length(degree), ncol = length(den))

  best <- list()

  for(d in 1:length(degree)){

    if(print_models){print(paste("degree:", degree[d]))}


    for(k in 1:length(den)){

      if(print_models){print(paste("knot position:  t_p /", den[k]))}

      n_knots   <- ifelse(!is.na(den[k]), 1,0)
      npar_spln <- degree[d] + n_knots + 1
      npar      <- npar_spln + 2


      loglikelihood <-  function(param){

        alpha      <- param[1]
        knot       <- param[2]
        gamma      <- param[3:(3+npar_spln-1)]


        lik <- 0
        for (j in 1:length(data)) {

          kk <- data[[j]][,2]

          tt <- log(t)[t<knot]
          rr <- log(t)[t>=knot]


          if(n_knots>0){
            mat <- bSpline(tt, intercept = TRUE,  degree = degree[d], knots = log(knot/den[k]))
          }else{
            mat <- bSpline(tt, intercept = TRUE,  degree = degree[d])
          }


          mat <- rbind(mat, matrix(0, nrow= length(rr)  ,ncol = ncol(mat)))
          mat <- cbind(mat, c(rep(0,length(tt)), rep(1,length(rr))) )


          if(length(kk)>=knot){
            kk1  <- sum(kk[knot:length(kk)])
            l1   <- length(knot:length(kk))
            indicator <- 1
            mat2 <- mat[1:(knot-1),]
            kk2  <- kk[1:(knot-1)]

          }else{
            kk1  <- 0
            l1   <- 0
            indicator <- 0
            mat2 <- mat[1:length(kk),]
            kk2  <- kk
          }


          lik <- lik + indicator*(lgamma(alpha+kk1)-lgamma(alpha)  +
                                    alpha*log(alpha) + kk1*(gamma[length(gamma)] )-
                                    (alpha+kk1)*log(alpha+l1*exp(gamma[length(gamma)])))+

            sum(
              lgamma(alpha+kk2)-lgamma(alpha)+
                alpha*log(alpha) +
                kk2*(mat2%*%c(gamma,gamma[length(gamma)]) )-
                (alpha+kk2)*log(alpha+exp(mat2%*%c(gamma,gamma[length(gamma)])))



            )



        }
        return(-lik)
      }




      par  <- matrix(NA, nrow = length(tp_start), ncol = npar)
      conv <- c()
      lik <- c()



      for(w in 1:length(tp_start)){



        op4 <- nloptr::bobyqa(x0 = c(1,tp_start[w], rep(1,npar-2)),fn =  loglikelihood,
                              control=list(maxeval=max_eval, xtol_rel=tol_rel),
                              lower =c(0.001, low_tp , rep(-Inf, npar-2)),
                              upper =c(Inf , t_int, rep(Inf, npar-2)))


        lik[w]  <- op4$value
        par[w,] <- op4$par


      }


      par_sel_tmp <- par[which.min(lik),]
      bic     <- npar*log(length(data)) + 2*min(lik)

      results[[d]][[k]] <- c( bic, par_sel_tmp )
      names(results) <- paste("degree:", degree)

      bic_mat[d,k] <- bic



    }


  }


  ind <- which(bic_mat == min(bic_mat), arr.ind = TRUE)


  if(!is.na(den[ind[2]])){
    best[[1]] <- paste("degree:", degree[ind[1]], ", knot position: t_p /", den[ind[2]])
  }else{
    best[[1]] <- paste("degree:", degree[ind[1]], ", no internal knots")

  }



  best[[2]] <- t(matrix(results[[ind[1]]][[ind[2]]][-1]))
  best[[3]] <- c(degree[ind[1]], den[ind[2]])
  colnames(best[[2]]) <- c("alpha", "t_p", paste("spline", 1:(length(best[[2]])-2)) )

  names(best) <- c("Chosen model", "Parameters", "Degree and knot")

  return(list(results=results, best_model=best, data = data, t_int = t_int ))



}














#' Recruitment predictions via the tPG model
#'
#' @param est Output from the \code{tPG_est} function
#' @param tpred Future time point for predictions. It can be a single value or a grid of time points.
#' @param start_init Starting times of the centers already activated.
#' @param N Optional; target sample size. If specified, \code{tpred} must be a grid of values.
#' @param start_add Starting times of centers not yet activated.
#' @param conf_level Confidence level of the credible interval.
#'
#' @return This function returns multiple objects, including the (additional) predicted enrollments and associated credible intervals from the current time to the values
#' specified in \code{tpred}. If \code{N} is specified, the function also returns the estimated time left to recruit the target sample size.
#' @export
#'
#' @examples \dontrun{sim_data <- sim_rec(alpha = 1, p1 = 10, p2 = 0.15, c2 = 0.3, c1 = 0.1,
#' type = "cdf", t_p = 150, Tmax = 200, s_time = 1:100, C_act = 30, random_init = TRUE)
#'
#' est <- tPG_est(data = sim_data$data, t_int = 200 , tp_start = c(20, 60, 100, 140), degree = c(2, 3), den = c(NA, 2, 3) )
#' pred <- tPG_pred(est = est, start_init = sim_data$start, N=5000, tpred = 201:600, conf_level = 0.95)}
tPG_pred <- function(est ,  tpred, start_init, N = NULL, start_add = NULL, conf_level = 0.95){

  data <-  est$data
  t_int <- est$t_int

  start   <- c(start_init, start_add)
  C_tot   <- length(start)
  C_act   <- length(start_init)
  par_sel <- est$best_model$Parameters
  t <- 1:t_int

  npar    <- length(par_sel)
  spln    <- par_sel[2]
  cutoff  <- ceiling(spln)

  if(cutoff>t_int){cutoff=t_int}

  index <- (t_int-start+1)-cutoff>0


  nj_t_int <- c()
  timespan <- c()
  for(j in 1:sum(index)){

    nj_t_int[j] <- sum(data[[j]][,2][cutoff:(length(data[[j]][,2]))])
    timespan[j] <- length(cutoff:(t_int-start[j]))
  }

  alpha_m <- par_sel[1]

  tt <- log(t)[t<spln]
  rr <- log(t)[t>=spln]


  if (!is.na(est$best_model$`Degree and knot`[2])){
    mat <- bSpline(tt, intercept = TRUE, knots = log(spln/est$best_model$`Degree and knot`[2]), degree =est$best_model$`Degree and knot`[1])
  }else{mat <- bSpline(tt, intercept = TRUE, degree = est$best_model$`Degree and knot`[1])}


  mat <- rbind(mat, matrix(0, nrow= length(rr)  ,ncol = ncol(mat)))
  mat <- cbind(mat, c(rep(0,length(tt)), rep(1,length(rr))) )

  lin <- mat%*%c(par_sel[3:(npar)], par_sel[(npar)])

  beta_m  <- alpha_m/exp(lin)[cutoff]

  a       <-  (alpha_m + nj_t_int)
  b       <-  (timespan + beta_m)


  beta_m2 <- alpha_m/exp(lin)

  sum_p <- sum(index==TRUE)

  tau   <- tpred-t_int



  m <- c()
  v <- c()
  ci <- matrix(NA, ncol=2, nrow = length(tau))

  for(k in 1:length(tau)){


    if(sum(index==TRUE)==C_tot){
      m2_2  <- 0
      m2v_2 <- 0
      timespan_2 <- 0
    }else{

      m2_2 <- c()
      m2v_2 <- c()
      timespan_2 <- c()

      for(j in (sum(index==TRUE)+1):C_tot){


        if(t_int-start[j]+tau[k]>0 & (max(1, t_int-start[j]+1)<=(min(cutoff-1, t_int-start[j]+tau[k])))){

          m2_2[j-(sum(index==TRUE))]  <- sum(alpha_m / beta_m2[max(1, t_int-start[j]+1):(min(cutoff-1, t_int-start[j]+tau[k]))])

          m2v_2[j-(sum(index==TRUE))] <- sum(alpha_m / (beta_m2[max(1, t_int-start[j]+1):(min(cutoff-1, t_int-start[j]+tau[k]))]^2))

          if(cutoff<=(tpred[k]-start[j])){

            timespan_2[j-(sum(index==TRUE))] <- length(cutoff:(tpred[k]-start[j])   )

          }else{timespan_2[j-(sum(index==TRUE))]=0}

        }else{

          m2_2[j-(sum(index==TRUE))] <- 0
          timespan_2[j-(sum(index==TRUE))] <- 0
          m2v_2[j-(sum(index==TRUE))] <- 0


        }






      }

    }


    m[k] <- tau[k]*sum(a/b) + sum(m2_2) + sum(timespan_2*(alpha_m/beta_m))
    v[k]    <- tau[k]*sum(a/b)+ sum(m2_2) + sum(timespan_2*(alpha_m/beta_m))+
      (tau[k]^2)*sum(a/b^2)+sum(m2v_2) +  sum((timespan_2^2)*(alpha_m/(beta_m^2)))

    ci[k,] <- m[k] + c(qnorm((1-conf_level)/2), -qnorm((1-conf_level)/2))*sqrt(v[k])



  }


  out_1 <- data.frame(m = m, ci_l = ci[,1], ci_u = ci[,2], t=tpred )
  if(is.numeric(N)){

    count <- NULL
    for(i in 1:C_act){
      count <- c(count, data[[i]][,2])
    }

    N_rem <- N-sum(count)

    out_2 <- c(which(m      > N_rem)[1],
               which(ci[,2]> N_rem)[1],
               which(ci[,1]> N_rem)[1])

    return(list(n_pred = out_1 , T_pred = out_2, N = N, t_int = t_int, start_add = start_add, start_init = start_init,
                data = data) )

  }else(  return(list(n_pred = out_1 , N = N, t_int = t_int, start_add = start_add, start_init = start_init,
                      data = data) )
  )



}











#' Display forecasted enrollments
#'
#' @param pred Output from the \code{tPG_pred} function.
#' @param xlim Standard \code{plot} argument
#' @param ylim Standard \code{plot} argument
#' @param print_T Logical. If TRUE, point estimate and CrIs for the time to achieve the sample size are printed.
#' @param main Standard \code{plot} argument
#' @param xlab Standard \code{plot} argument
#' @param ylab Standard \code{plot} argument
#'
#' @return A plot of the recruitment process.
#' @export
#'
#' @examples sim_data <- sim_rec(alpha = 1, p1 = 10, p2 = 0.15, c2 = 0.3, c1 = 0.1,
#' type = "cdf", t_p = 150, Tmax = 200, s_time = 1:100, C_act = 30, random_init = TRUE)
#'
#' est <- tPG_est(data = sim_data$data, t_int = 200 , tp_start = c(20, 60, 100, 140), degree = c(2, 3), den = c(NA, 2, 3) )
#' pred <- tPG_pred(est = est, start_init = sim_data$start, N=5000, tpred = 201:900, conf_level = 0.95)
#' tPG_plot(pred = pred, xlim=c(0,550), ylim=c(0,5500), print_T = TRUE, xlab = "Days", ylab = "Enrollments", main = "Forecasted recruitments")
tPG_plot <- function(pred,  xlim, ylim, print_T = F, main = NULL, xlab = NULL, ylab = NULL){

  start_init <- pred$start_init
  data <- pred$data
  C_act   <- length(start_init)
  count <- NULL
  time  <- NULL
  for(i in 1:C_act){


    count <- c(count, data[[i]][,2])

    time <-  c(time, data[[i]][,1] + (start_init[i]-1))

  }


  rec_mat <- cbind(time, count)
  rec_mat <- rec_mat[order(rec_mat[,1]),]



  plot(rec_mat[,1], cumsum(rec_mat[,2]), type="l", ylim = ylim, xlim=xlim, main = main,
       xlab = xlab, ylab = ylab )
  lines(pred$n_pred$t, pred$n_pred$m + sum(rec_mat[,2]) , col= "red"   )
  lines(pred$n_pred$t, pred$n_pred$ci_l + sum(rec_mat[,2]) , col= "red" , lwd="1", lty=2  )
  lines(pred$n_pred$t, pred$n_pred$ci_u + sum(rec_mat[,2]) , col= "red"  , lwd="1", lty=2 )

  points(c(start_init, pred$start_add), rep(0,length(c(start_init, pred$start_add))), col = "red", pch = 3)


  if(print_T){

    lines(-10000:pred$T_pred[3]+pred$t_int, rep(pred$N, length(-10000:pred$T_pred[3]+pred$t_int)) ,   lwd="1", lty=2  )

    lines(rep(pred$T_pred[1]+pred$t_int, length(0:pred$N)), 0:pred$N ,   lwd="1", lty=2  )
    lines(rep(pred$T_pred[3]+pred$t_int, length(0:pred$N)), 0:pred$N ,  lwd="1", lty=2  )
    lines(rep(pred$T_pred[2]+pred$t_int, length(0:pred$N)), 0:pred$N ,  lwd="1", lty=2  )

  }


}














#' Estimating function of the tPG model parameters with covariates
#'
#' @param data Recruitment data of the study formatted as follows: a list where each element contains the recruitment data of one center formatted as a data frame with two columns;
#' the first column is the grid of recruitment days (or another time unit) for the concerned center, and the second column indicates the respective number of enrollments on that day.
#' @param X Matrix of center-specific covariates where each row is a covariate.
#' The centers represented by row must be in the same chronologial order used in \code{data}.
#' @param t_int Interim time of the analysis with respect to the start of the first center.
#' @param tp_start Starting values for the plateau point to be used in the optimization function. The likelihood is optimized for all these starting values and the the estimated parameters
#' corresponding to the highest likelihood are used.
#' @param degree Candidate degree values for the polynomial.
#' @param den Candidate internal knot placement, expressed as a function of the estimated plateau point, i.e., t_p/\code{den}
#' @param max_eval Maximum number of iterations of the optimization function. Default is 100000.
#' @param tol_rel Tolerance of the optimization function. Default is 1.e-8.
#' @param print_models Logical. Print iterations across internal knots and polynomial degrees candidates.
#' @param low_tp Lower bound for the plateau point to be used in the optimization function.
#'
#' @return This function returns multiple objects which include the estimates of the parameters under each model and the best-fitting model selected via BIC.
#' @export
#'@example
#' sim_data <- sim_rec(alpha = 1, p1 = 10, p2 = 0.15, c2 = 0.3, c1 = 0.1,
#'   type = "cdf", t_p = 150, Tmax = 200, s_time = 1:100, C_act = 30, random_init = TRUE) # Note that there are no covariates in the data generating function
#' X <- matrix(rbinom(60, size = 1, prob=0.5), ncol=2)
#' est <- tPG_est_cov(data = sim_data$data, X=X, t_int = 200 , tp_start = c(20, 60, 100, 140), degree = c(2, 3), den = c(NA, 2) )
#'
tPG_est_cov <- function(data, X = X, t_int, tp_start, degree, den, max_eval=100000, tol_rel=1.e-8, print_models = T , low_tp = 6){




  t <- 1:t_int


  results <- rep(list(list()),length(degree))
  bic_mat <- matrix(NA, nrow = length(degree), ncol = length(den))

  best <- list()

  for(d in 1:length(degree)){

    if(print_models){print(paste("degree:", degree[d]))}


    for(k in 1:length(den)){

      if(print_models){print(paste("knot position:  t_p /", den[k]))}

      n_knots   <- ifelse(!is.na(den[k]), 1,0)
      npar_spln <- degree[d] + n_knots + 1
      npar      <- npar_spln + 2 + ncol(X)


      loglikelihood <-  function(param){

        alpha      <- param[1]
        knot       <- param[2]
        gamma      <- param[3:(3+npar_spln-1)]
        x_coef     <- param[(3+npar_spln):npar]

        lik <- 0
        for (j in 1:length(data)) {

          kk <- data[[j]][,2]

          tt <- log(t)[t<knot]
          rr <- log(t)[t>=knot]


          if(n_knots>0){
            mat <- bSpline(tt, intercept = T,  degree = degree[d], knots = log(knot/den[k]))
          }else{
            mat <- bSpline(tt, intercept = T,  degree = degree[d])
          }


          mat <- rbind(mat, matrix(0, nrow= length(rr)  ,ncol = ncol(mat)))
          mat <- cbind(mat, c(rep(0,length(tt)), rep(1,length(rr))) )


          if(length(kk)>=knot){
            kk1  <- sum(kk[knot:length(kk)])
            l1   <- length(knot:length(kk))
            indicator <- 1
            mat2 <- mat[1:(knot-1),]
            kk2  <- kk[1:(knot-1)]

          }else{
            kk1  <- 0
            l1   <- 0
            indicator <- 0
            mat2 <- mat[1:length(kk),]
            kk2  <- kk
          }


          lik <- lik + indicator*(lgamma(alpha+kk1)-lgamma(alpha)  +
                                    alpha*log(alpha) + kk1*(gamma[length(gamma)] + sum(c(X[j,])*c(x_coef)) )-
                                    (alpha+kk1)*log(alpha+l1*exp(gamma[length(gamma)] + sum(c(X[j,])*c(x_coef)))))+

            sum(
              lgamma(alpha+kk2)-lgamma(alpha)+
                alpha*log(alpha) +
                kk2*(mat2%*%c(gamma,gamma[length(gamma)]) + sum(c(X[j,])*c(x_coef)))-
                (alpha+kk2)*log(alpha+exp(mat2%*%c(gamma,gamma[length(gamma)])+ sum(c(X[j,])*c(x_coef))))



            )



        }
        return(-lik)
      }




      par  <- matrix(NA, nrow = length(tp_start), ncol = npar)
      conv <- c()
      lik <- c()



      for(w in 1:length(tp_start)){



        op4 <- nloptr::bobyqa(x0 = c(1,tp_start[w], rep(1,npar-2)),fn =  loglikelihood,
                              control=list(maxeval=max_eval, xtol_rel=tol_rel),
                              lower =c(0.001, low_tp , rep(-Inf, npar-2)),
                              upper =c(Inf , t_int, rep(Inf, npar-2)))


        lik[w]  <- op4$value
        par[w,] <- op4$par


      }


      par_sel_tmp <- par[which.min(lik),]
      bic     <- npar*log(length(data)) + 2*min(lik)

      results[[d]][[k]] <- c( bic, par_sel_tmp )
      names(results) <- paste("degree:", degree)

      bic_mat[d,k] <- bic



    }


  }


  ind <- which(bic_mat == min(bic_mat), arr.ind = TRUE)


  if(!is.na(den[ind[2]])){
    best[[1]] <- paste("degree:", degree[ind[1]], ", knot position: t_p /", den[ind[2]])
  }else{
    best[[1]] <- paste("degree:", degree[ind[1]], ", no internal knots")

  }



  best[[2]] <- t(matrix(results[[ind[1]]][[ind[2]]][-1]))
  best[[3]] <- c(degree[ind[1]], den[ind[2]])
  colnames(best[[2]]) <- c("alpha", "t_p", paste("spline", 1:(length(best[[2]])-2-ncol(X))), paste("X_coef", 1:(ncol(X)) ))

  names(best) <- c("Chosen model", "Parameters", "Degree and knot")

  return(list(results=results, best_model=best, data = data, t_int = t_int ))



}



#' Recruitment predictions via the tPG model with covariates
#'
#' @param est Output from the \code{tPG_est} function
#' @param X Matrix of center-specific covariates where each row is a covariate.
#' The centers represented by row must be in the same chronologial order used in \code{data}.
#' @param tpred Future time point for predictions. It can be a single value or a grid of time points.
#' @param start_init Starting times of the centers already activated.
#' @param N Optional; target sample size. If specified, \code{tpred} must be a grid of values.
#' @param start_add Starting times of centers not yet activated.
#' @param conf_level Confidence level of the credible interval.
#'
#' @return This function returns multiple objects, including the (additional) predicted enrollments and associated credible intervals from the current time to the values
#' specified in \code{tpred}. If \code{N} is specified, the function also returns the estimated time left to recruit the target sample size.
#' @export
#'
#' @example
#' sim_data <- sim_rec(alpha = 1, p1 = 10, p2 = 0.15, c2 = 0.3, c1 = 0.1,
#'   type = "cdf", t_p = 150, Tmax = 200, s_time = 1:100, C_act = 30, random_init = TRUE) # Note that there are no covariates in the data generating function
#' X <- matrix(rbinom(60, size = 1, prob=0.5), ncol=2)
#' est <- tPG_est_cov(data = sim_data$data, X=X, t_int = 200 , tp_start = c(20, 60, 100, 140), degree = c(2, 3), den = c(NA, 2) )
#'
#' pred <- tPG_pred_cov(est = est, X=X,start_init = sim_data$start, N=5000, tpred = 201:600, conf_level = 0.95)
#'
tPG_pred_cov <- function(est , X = X, tpred, start_init, N = NULL, start_add = NULL, conf_level = 0.95){

  data <-  est$data
  t_int <- est$t_int

  start   <- c(start_init, start_add)
  C_tot   <- length(start)
  C_act   <- length(data)
  par_sel <- est$best_model$Parameters
  t <- 1:t_int

  npar    <- length(par_sel)
  spln    <- par_sel[2]
  cutoff  <- ceiling(spln)

  if(cutoff>t_int){cutoff=t_int}

  index <- (t_int-start+1)-cutoff>0


  nj_t_int <- c()
  timespan <- c()
  for(j in 1:sum(index)){

    nj_t_int[j] <- sum(data[[j]][,2][cutoff:(length(data[[j]][,2]))])
    timespan[j] <- length(cutoff:(t_int-start[j]))
  }


  alpha_m <- par_sel[1]

  tt <- log(t)[t<spln]
  rr <- log(t)[t>=spln]


  if (!is.na(est$best_model$`Degree and knot`[2])){
    mat <- bSpline(tt, intercept = T, knots = log(spln/est$best_model$`Degree and knot`[2]), degree =est$best_model$`Degree and knot`[1])
  }else{mat <- bSpline(tt, intercept = T, degree = est$best_model$`Degree and knot`[1])}


  mat <- rbind(mat, matrix(0, nrow= length(rr)  ,ncol = ncol(mat)))
  mat <- cbind(mat, c(rep(0,length(tt)), rep(1,length(rr))) )



  lin <- matrix(NA, ncol = C_tot, nrow = nrow(mat))

  for(j in 1:C_tot){

    lin[,j] <- mat%*%c(par_sel[3:(npar-ncol(X))], par_sel[(npar-ncol(X))])+
      sum(X[j,]*par_sel[(npar-ncol(X)+1):npar])

  }

  beta_m  <- alpha_m/exp(lin)[cutoff,1:sum(index==T)]


  a       <-  (alpha_m + nj_t_int)
  b       <-  (timespan + beta_m)



  sum_p <- sum(index==T)

  tau   <- tpred-t_int



  m <- c()
  v <- c()
  ci <- matrix(NA, ncol=2, nrow = length(tau))

  for(k in 1:length(tau)){


    if(sum(index==T)==C_tot){
      m2_2  <- 0
      m2v_2 <- 0
      timespan_2 <- 0
      beta_m3  <- 0

      m_a <- tau*sum(a/b) + sum(m2_2)
      v    <- tau*sum(a/b)+ sum(m2_2) +   (tau^2)*sum(a/b^2)+sum(m2v_2)
    }else{

      m2_2 <- c()
      m2v_2 <- c()
      timespan_2 <- c()

      for(j in (sum(index==T)+1):C_tot){

        beta_m2 <- (alpha_m/exp(lin))[,j]

        if(t_int-start[j]+tau[k]>0  &  (max(1, t_int-start[j]+1)<=(min(cutoff-1, t_int-start[j]+tau[k])))   ){


          m2_2[j-(sum(index==T))]  <- sum(alpha_m / beta_m2[max(1, t_int-start[j]+1):(min(cutoff-1, t_int-start[j]+tau[k]))])

          m2v_2[j-(sum(index==T))] <- sum(alpha_m / (beta_m2[max(1, t_int-start[j]+1):(min(cutoff-1, t_int-start[j]+tau[k]))]^2))

          if(cutoff<=(tpred[k]-start[j])){

            timespan_2[j-(sum(index==T))] <- length(cutoff:(tpred[k]-start[j])   )

          }else{timespan_2[j-(sum(index==T))]=0}

        }else{

          m2_2[j-(sum(index==T))] <- 0
          timespan_2[j-(sum(index==T))] <- 0
          m2v_2[j-(sum(index==T))] <- 0

        }






      }

    }

    beta_m3  <- alpha_m/exp(lin)[cutoff,ifelse(sum(index==T)==C_tot, 1,(sum(index==T)+1)):C_tot]

    m[k] <- tau[k]*sum(a/b) + sum(m2_2) + sum(timespan_2*(alpha_m/beta_m3))
    v[k]    <- tau[k]*sum(a/b)+ sum(m2_2) + sum(timespan_2*(alpha_m/beta_m3))+
      (tau[k]^2)*sum(a/b^2)+sum(m2v_2) +  sum((timespan_2^2)*(alpha_m/(beta_m3^2)))

    ci[k,] <- m[k] + c(qnorm((1-conf_level)/2), -qnorm((1-conf_level)/2))*sqrt(v[k])



  }


  out_1 <- data.frame(m = m, ci_l = ci[,1], ci_u = ci[,2], t=tpred )
  if(is.numeric(N)){

    count <- NULL
    for(i in 1:C_act){
      count <- c(count, data[[i]][,2])
    }

    N_rem <- N-sum(count)

    out_2 <- c(which(m      > N_rem)[1],
               which(ci[,2]> N_rem)[1],
               which(ci[,1]> N_rem)[1])

    return(list(n_pred = out_1 , T_pred = out_2, N = N, t_int = t_int, start_add = start_add, start_init = start_init,
                data = data) )

  }else(  return(list(n_pred = out_1 , N = N, t_int = t_int, start_add = start_add, start_init = start_init,
                      data = data) )
  )



}


















