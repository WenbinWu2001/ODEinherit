smoother_SS <- function(obs_time, Y, tt){
  # smoothing spline estimates (KernelODE step 1)
  # Y is (n,p) containing the observed time series.
  # tt is the (finer) time grid to evaluate the smoothing spline.
  p <- ncol(Y)
  times_e <- c(0,tt)  # also fit an initial value
  init_vals_smth <- rep(NA, p)
  yy_smth <- matrix(NA, length(tt), p)  # on tt
  deriv_smth <- matrix(NA, length(tt), p)  # on tt
  for(j in 1:p){
    SS_model <- stats::smooth.spline(obs_time, Y[,j], all.knots = T)  # a cubic smoothing spline where all points in `obs_time` are used as knots
    gcvscore <- SS_model$cv.crit
    pred_times_e <- stats::predict(SS_model, times_e)$y

    init_vals_smth[j] <- pred_times_e[1]
    yy_smth[,j] <- pred_times_e[-1]
    deriv_smth[,j] <- (stats::predict(SS_model,times_e,deriv=1)$y)[-1]
  }

  return (list(init_vals_smth = init_vals_smth,
               yy_smth = yy_smth,
               deriv_smth = deriv_smth))
}






# smoothX <- function(observations,
#                     times_e,
#                     deg,
#                     h=NULL,
#                     maxk=5000,
#                     type_smooth=c("local polynomial","smoothspline" ,"linear interpolation"),
#                     type_data="Replicates"){
#   # function in GRADE to smooth the data
#   # TODO: Add citation, add notes on input and output
#   # `observation` is a list of R elements, each being an n x (1+p) matrix where the first column being the time points for the observations.
#   Xprime<-Xhat<-list()
#   gcvscore<-list()
#   if(type_data=="Replicates"){
#     obs<-do.call(rbind,observations)
#     p<-dim(obs)[2]-1
#
#     Xprime[[1]]<-Xhat[[1]]<-matrix(0,nrow=length(times_e),ncol= p)
#
#     if( type_smooth=="linear interpolation") {
#       for(i in 1:p){Xhat[[1]][,i]<-approx(x=obs[,1],y=obs[,i+1] ,xout=times_e,rule=2)$y
#       gcvscore[[1]]<-NULL
#       Xprime[[1]]<-NULL
#       deg<-NULL
#       }
#     } else if(type_smooth=="local polynomial"){
#       gcvscore[[1]]<-matrix(0,nrow=length(h),ncol= p)
#
#       # Get the LP estimates
#       for(i in 1:p){
#         # For the lp estimates, we pick the bandwidth that minimize the generalized cross-validation
#         # for each node.
#         for(j in 1:length(h)){
#           h_temp<-h[j]
#           temp <- locfit(obs[,i+1]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)
#           gcvscore[[1]][j,i]<-gcv(obs[,i+1]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)[4]
#         }
#         h_temp<-h[which.min(gcvscore[[1]][,i])]
#         temp <- locfit(obs[,i]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)
#         Xhat[[1]][,i]<-stats::predict(temp,newdata=times_e)
#         Xprime[[1]][,i]<-  NULL
#       }
#
#     } else if (type_smooth == "smoothspline"){
#       for(i in 1:p){
#         gcvscore[[1]]<-numeric(p)
#
#         temp <- smooth.spline(obs[,1], obs[,i+1], all.knots = T)
#         gcvscore[[1]][i]<-temp$cv.crit
#         Xhat[[1]][,i]<-stats::predict(temp,times_e)$y
#         Xprime[[1]][,i]<-  stats::predict(temp,times_e,deriv=1)$y
#       }
#     }
#   } else if (type_data=="Perturbations"){
#     # Similar to the previous function, but create a list for every experiment
#     R<-length(observations)
#     for(r in 1:R){  # smooth for each experiment separately
#       obs<-observations[[r]]
#       p<-dim(obs)[2]-1
#
#       Xprime[[r]]<-Xhat[[r]]<-matrix(0,nrow=length(times_e),ncol= p)
#
#       if( type_smooth=="linear interpolation") {
#         for(i in 1:p){
#           Xhat[[r]][,i]<-approx(x=obs[,1],y=obs[,i+1] ,xout=times_e,rule=2)$y
#           gcvscore[[r]]<-NULL
#           Xprime[[r]]<-NULL
#           deg<-NULL
#         }
#       } else if(type_smooth=="local polynomial"){
#         # Get the LP estimates
#         gcvscore[[r]]<-matrix(0,nrow=length(h),ncol= p)
#         for(i in 1:p){
#           for(j in 1:length(h)){
#             h_temp<-h[j]
#             temp <- locfit(obs[,i+1]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)
#             gcvscore[[r]][j,i]<-gcv(obs[,i+1]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)[4]
#           }
#           h_temp<-h[which.min(gcvscore[[r]][,i])]
#           temp <- locfit(obs[,i]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)
#           Xhat[[r]][,i]<-stats::predict(temp,newdata=times_e)
#           Xprime[[r]][,i]<-  NULL
#         }
#       } else if (type_smooth == "smoothspline"){
#         gcvscore[[r]]<-numeric(p)
#         for(i in 1:p){
#           temp <- smooth.spline(obs[,1], obs[,i+1], all.knots = T)
#           gcvscore[[r]][i]<-temp$cv.crit
#           Xhat[[r]][,i]<-stats::predict(temp,times_e)$y
#           Xprime[[r]][,i]<-  stats::predict(temp,times_e,deriv=1)$y
#         }
#       }
#     }
#   }
#   return(list(Xhat=Xhat,Xprime=Xprime,times_e=times_e,gcvscore=gcvscore,type_smooth=type_smooth))
# }


# .SS_estimation <- function(obs_time,
#                            Y,
#                            sp = 1e-4,
#                            SS_dof = 20,
#                            verbose = 0){
#   ## smoothing spline estimation
#   ## Y is centered and scaled data matrix or data.frame with shape nxp.
#   ## sp and SS_dof are arguments of s() that defines the smooth terms. sp is the smoothing parameter, SS_dof is the k argument (the dimension of the basis used to represent the smooth term).
#   ## Returns a list of mgcv::gam models.
#   if (verbose > 0) {cat("-------- estimating smoothing splines --------\n")}
#
#   SS_fit_list <- list()
#   df <- data.frame(cbind(Y, obs_time))
#   for (j in 1:ncol(Y)){
#     # SS_fit_list[[j]] <- npregfast::frfast(as.formula(paste0(colnames(df)[j], "~s(obs_time, m = 2, k =3)")), data = df, nboot = 100, smooth = "splines", cluster = TRUE, ncores = 2)
#     # SS_fit_list[[j]] <- mgcv::gam(as.formula(paste0(colnames(df)[j], "~s(obs_time, sp = ", sp ,", k = ", SS_dof, ", bs = 'tp')")), data = df)
#     SS_fit_list[[j]] <- mgcv::gam(as.formula(paste0(colnames(df)[j], "~s(obs_time, sp = 1e-4, k = 20, bs = 'tp')")), data = df)
#   }
#   return (SS_fit_list)
# }
#
#
# .SS_interpolation <- function(SS_fit_list, tt){
#   yy_smth <- do.call(cbind, lapply(SS_fit_list, function(SS_fit){stats::predict(SS_fit, newdata = data.frame(obs_time = tt))}))
#   # yy_smth <- do.call(cbind, lapply(SS_fit_list, function(SS_fit){stats::predict(SS_fit, newdata = data.frame(obs_time = tt))$Estimation[,'Pred']}))
#   return (yy_smth)
# }
