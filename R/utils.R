prepare <- function(x) {
  # This function is used for prepare standata for normal and lognormal model
  tmp_data <- na.omit(x)  |>
    filter(rt > 0)

  n_subj <- length(unique(tmp_data$subj_num))# Number of subjects

  n_cond <- 2 # Number of conditions

  n_time <- 2 # Number of timepoints

  # Determine number of trials within subjects, conditions, and timepoints
  T_subj_all <- tmp_data |>
    group_by(subj_num, Condition, test_index)  |>
    summarize(n_trials = n()) # T_subj_all[subj_num, Condition, test_index, trial数]

   T_max <-  max(T_subj_all$n_trials)

  # Create RT data array for stan; dims = (subject, condition, time, trial)
  RT <- correct <- array(0, dim = c(n_subj, n_cond, n_time, T_max))
  T_subj <- RT_min <- array(NA, dim = c(n_subj, n_cond, n_time)) # number of trials, min RT
  for (i in 1:n_subj) {
    # Number of trials per condition/timpoint
    c0_t1 <- T_subj[i, 1, 1] <- with(T_subj_all, n_trials[subj_num==i & Condition==0 & test_index==1])
    c2_t1 <- T_subj[i, 2, 1] <- with(T_subj_all, n_trials[subj_num==i & Condition==2 & test_index==1])
    c0_t2 <- T_subj[i, 1, 2] <- with(T_subj_all, n_trials[subj_num==i & Condition==0 & test_index==2])
    c2_t2 <- T_subj[i, 2, 2] <- with(T_subj_all, n_trials[subj_num==i & Condition==2 & test_index==2])

    # RTs and correctness for congruent condition at time 1
    RT[i, 1, 1, 1:c0_t1]      <- with(tmp_data, rt[subj_num==i & Condition==0 & test_index==1])
    RT_min[i, 1, 1]           <- min(RT[i, 1, 1, 1:c0_t1])
    correct[i, 1, 1, 1:c0_t1] <- with(tmp_data, acc[subj_num==i & Condition==0 & test_index==1])
    # Choice and RTs for incongruent condition at time 1
    RT[i, 2, 1, 1:c2_t1]      <- with(tmp_data, rt[subj_num==i & Condition==2 & test_index==1])
    RT_min[i, 2, 1]           <- min(RT[i, 2, 1, 1:c2_t1])
    correct[i, 2, 1, 1:c2_t1] <- with(tmp_data, acc[subj_num==i & Condition==2 & test_index==1])
    # Choice and RTs for congruent condition at time 2
    RT[i, 1, 2, 1:c0_t2]      <- with(tmp_data, rt[subj_num==i & Condition==0 & test_index==2])
    RT_min[i, 1, 2]           <- min(RT[i, 1, 2, 1:c0_t2])
    correct[i, 1, 2, 1:c0_t2] <- with(tmp_data, acc[subj_num==i & Condition==0 & test_index==2])
    # Choice and RTs for incongruent condition at time 2
    RT[i, 2, 2, 1:c2_t2]      <- with(tmp_data, rt[subj_num==i & Condition==2 & test_index==2])
    RT_min[i, 2, 2]           <- min(RT[i, 2, 2, 1:c2_t2])
    correct[i, 2, 2, 1:c2_t2] <- with(tmp_data, acc[subj_num==i & Condition==2 & test_index==2])
  }

  # Stan-ready data list
  list(N       = n_subj,
       N_cond  = n_cond,
       N_time  = n_time,
       T_max   = T_max,
       T_subj  = T_subj,
       RT      = RT,
       RT_min  = RT_min,
       correct = correct)
}


prepare_wiener <- function(x){
  # This function is used for prepare standata for Wiener model

  tmp_data <- na.omit(x)  |>
    filter(rt > 0)

  tmp_data <- tmp_data |> mutate(
    subj_num = as.integer(factor(user_id, levels = unique(user_id)))
  )

  n_subj <- length(unique(tmp_data$subj_num))#amount of subjects

  n_cond <- 2 # Number of conditions

  # n_time <- 2 # Number of timepoints

  # Determine number of trials within subjects, conditions, and timepoints
  T_subj_all <- tmp_data |>
    group_by(subj_num, Condition,acc)  |>
    summarize(n_trials = n()) # T_subj_all[subj_num, Condition, test_index, trial数]

  T_max <-  max(T_subj_all$n_trials)

  #SUBJECT_ID, REACTION TIME, STIMULI SIDE (alternatively ANSWER SIDE), CORRECTNESS OF RESPONSE
  Nu <- tmp_data |> filter(acc==1) |> group_by(subj_num) |> summarize(n_trials = n()) |> pull(n_trials) #amount of upper boundary responses per subject
  Nl <- tmp_data |> filter(acc==0) |> group_by(subj_num) |> summarize(n_trials = n()) |> pull(n_trials) #amount of lower boundary responses per subject
  T_max_u <- max(Nu)
  T_max_l <- max(Nl)
  minRT <- tmp_data |> group_by(subj_num) |> summarize(value = min(rt)) |> pull(value) #min RT per subject

  # Create RT data array for stan; dims = (subject, condition, time, trial)
  # Reaction times for upper and lower boundary responses, PADDED matrices
  RTu <- array(0, dim = c(n_subj, n_cond, T_max_u)) # 构建一个n_subj x n_cond x T_max_u的0矩阵
  RTl <- array(0, dim = c(n_subj, n_cond, T_max_l)) #
  T_subj <- array(NA, dim = c(n_subj, n_cond)) # 构建一个n_subj x n_cond的矩阵

  # Store each subjects' reaction time data

  #ADAPT: DEFINE THE CODING OF THE MODEL:
  #STIMULI CODING:  tmp$rt[tmp$SIMULI_SIDE==1 or 0]
  #ACCURACY CODING:  tmp$rt[tmp$CORRECTNESS==1 or 0]
  for (i in 1:n_subj){
    # Number of trials per condition/timpoint
    c0_t1_u <- T_subj[i, 1] <- with(T_subj_all, n_trials[subj_num==i & Condition==0 & acc==1])
    c0_t1_l <- T_subj[i, 1] <- with(T_subj_all, n_trials[subj_num==i & Condition==0 & acc==0])
    c2_t1_u <- T_subj[i, 2] <- with(T_subj_all, n_trials[subj_num==i & Condition==2 & acc==1])
    c2_t1_l <- T_subj[i, 2] <- with(T_subj_all, n_trials[subj_num==i & Condition==2 & acc==0])

    # RTs for congruent condition at time 1
    RTu[i, 1, 1:c0_t1_u]      <- with(tmp_data, rt[subj_num==i & Condition==0 & acc==1])
    RTl[i, 1, 1:c0_t1_l]      <- with(tmp_data, rt[subj_num==i & Condition==0 & acc==0])
    # RTs for incongruent condition at time 1
    RTu[i, 2, 1:c2_t1_u]      <- with(tmp_data, rt[subj_num==i & Condition==2 & acc==1])
    RTl[i, 2, 1:c2_t1_l]      <- with(tmp_data, rt[subj_num==i & Condition==2 & acc==0])
}
  RTbound <- 0.1

  # List of data sent to Stan
  forstan <- list(
    N       = n_subj, # Number of subjects
    N_cond  = n_cond, # Number of conditions
    Nu_max  = T_max_u,  # Max (across subjects) number of upper boundary responses
    Nl_max  = T_max_l,  # Max (across subjects) number of lower boundary responses
    Nu      = Nu,       # Number of upper boundary responses for each subj
    Nl      = Nl,       # Number of lower boundary responses for each subj
    RTu     = RTu,      # upper boundary response times
    RTl     = RTl,      # lower boundary response times
    minRT   = minRT,    # minimum RT for each subject of the observed data
    RTbound = RTbound   # lower bound or RT across all subjects (e.g., 0.1 second)
  )

  return(list(forstan = forstan))
}

organize_samples <- function(samples) {
  # This function is used for reorganize stan posterior samples, for further calculation of rhat, ess, etc.
  # It will return a list contains different parameters' samples, each element is a matrix with dimension of (num_chains, num_iterations)
  num_iterations <- dim(samples)[1]
  num_chains <- dim(samples)[2]
  num_parameters <- dim(samples)[3]

  organized_list <- vector("list", length = num_parameters)

  for (param in 1:num_parameters) {
    matrix_for_param <- matrix(nrow = num_chains, ncol = num_iterations)

    for (chain in 1:num_chains) {
      matrix_for_param[chain, ] <- samples[, chain, param]
    }

    organized_list[[param]] <- matrix_for_param
  }

  names(organized_list) <- attr(samples, "dimnames")$parameters
  return(organized_list)
}


compute_rhat <- function(samples) {
  # Number of chains and samples per chain
  m <- ncol(samples)
  n <- nrow(samples)

  # Within chain variance
  W <- mean(apply(samples, 2, var))

  # Between chain variance
  B <- n * var(colMeans(samples))

  # Variance of the means
  var_hat <- (1 - 1/n) * W + 1/n * B

  # R-hat
  Rhat <- sqrt(var_hat / W)
  return(Rhat)
}


compute_ess <- function(samples) {
  # Number of chains and samples per chain
  m <- nrow(samples)  # Number of chains
  n <- ncol(samples)  # Number of samples per chain

  # Within chain variance
  W <- mean(apply(samples, 1, var))  # Apply variance along rows (chains)

  # Between chain variance
  B <- n * var(rowMeans(samples))  # Mean of each chain

  # Variance of the means
  var_hat <- (1 - 1/n) * W + 1/n * B

  # Effective Sample Size
  rho_hat_s <- sapply(1:(n-2), function(s) {
    mean(sapply(1:m, function(chain) {
      cor(samples[chain, 1:(n-s)], samples[chain, (1+s):n])
    }))
  })

  # Autocorrelation estimate
  rho_hat_even = 1 - (W / var_hat) * (1 + sum(rho_hat_s))

  # ESS
  ESS <- m * n / rho_hat_even
  return(ESS)
}





