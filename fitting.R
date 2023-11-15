library(tidyverse)
library(ggplot2)
library(rstan)
rm(list=ls())
requireNamespace("bit64")
# task-specific preprocessing functions
setwd("D:/Project/CAMP/behav-models")
source("./utils.R")
targets::tar_load(raw_data_parsed_224379118576069)

# prepare the data ----
data <- raw_data_parsed_224379118576069 |> #选游戏
  group_by(user_id) |> filter(n()==2) |> ungroup() |> # 选玩过两次以上的玩家
  arrange(user_id, game_time)  |> group_by(user_id)  |> mutate(test_index = row_number()) |> ungroup() |> # 选出每个玩家的第一次和第二次
  select(c("user_id","test_index","raw_parsed")) |> unnest(raw_parsed)|>  # 展开
  filter(all(device != "mouse"), .by = user_id) |>  # 选出所有没有用鼠标的玩家
  select(c("user_id","test_index","type","acc","rt")) |>
  filter(rt>=100 & rt<=2000) |> # 选出反应时在100-2000ms之间的数据, 62 trials <100, 42 trials>2000
  # maybe better use rt<=1500 according to Ratcliff & Mackoon (2007)
  mutate(rt = rt/1000)  # 转换为秒

data <- data |> mutate(
  subj_num = as.integer(factor(user_id, levels = unique(user_id))),
  Condition = ifelse(type == "congruent", 0, 2),
)


# partial_data <- data |> group_by(subj_num) |> filter(n()==400) |> ungroup() |>
#   filter(subj_num %in% c(2:11)) |>  mutate(subj_num = subj_num-1)

# Fit normal and lognormal model----
stan_data <- prepare(data)
saveRDS(stan_data,"./Data/Stan_ready_data.rds")

# seed for shifted_log_normal
tmp_seed  <- 1
tmp_adapt <- .9

tmp_seed  <- 43201
tmp_adapt <- .8

# fit normal model
fit <- stan(file ="normal_model.stan",
            data   = stan_data,
            iter   = 3000,
            warmup = 1000,
            thin   = 5,
            chains = 4,
            cores  = 4,
            seed   = tmp_seed,
            control = list(adapt_delta = tmp_adapt)
  )
saveRDS(fit, "./Results/fit_normal.rds")

# fit lognormal model
fit <- stan(file ="./stanModel/lognormal_model.stan",
            data   = stan_data,
            iter   = 3000,
            warmup = 1000,
     #       thin   = 5,
            chains = 4,
            cores  = 4,
            seed   = tmp_seed,
            control = list(adapt_delta = tmp_adapt)
)
saveRDS(fit, "./Results/fit_lognormal.rds")
# saveRDS(fit, "./Results/fit_lognormal_partial.rds")

# fit shifted-lognormal model
fit <- stan(file ="./stanModel/shifted_lognormal_model.stan",
            data   = stan_data,
            iter   = 3000,
            warmup = 1000,
            #       thin   = 5,
            chains = 4,
            cores  = 4,
            seed   = tmp_seed,
            control = list(adapt_delta = tmp_adapt)
)
saveRDS(fit, "./Results/fit_shifted_lognormal.rds")

# Fit wiener model ----
new_data <- data |> group_by(user_id, Condition)  |>
  summarize(mean_accuracy = mean(acc)) |>
  group_by(user_id) |>
  filter(all(mean_accuracy < 1)) |>
  ungroup() |>   select(user_id)  |>   distinct() |>
  inner_join(data,by="user_id")
stan_data <- prepare_wiener(new_data)
saveRDS(stan_data,"./Data/Stan_wiener_data.rds")

fit <- stan(file ="./stanModel/Wiener_model.stan",
            data   = stan_data,
            iter   = 3000,
            warmup = 1000,
            thin   = 5,
            chains = 4,
            cores  = 4,
            seed   = tmp_seed,
            control = list(adapt_delta = tmp_adapt)
)
saveRDS(fit, "./Results/fit_wiener.rds")
# saveRDS(fit, "./Results/fit_wiener_partial.rds")
