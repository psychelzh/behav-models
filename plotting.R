library(tidyverse)
library(ggplot2)
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

rm("raw_data_parsed_224379118576069")

data <- data |> mutate(
  subj_num = as.integer(factor(user_id, levels = unique(user_id))),
  Condition = ifelse(type == "congruent", 0, 2),
)


#

data |> filter(test_index==1,subj_num==1) |>  ggplot()+
  geom_histogram(aes(rt))+facet_wrap(~Condition)+
  xlabel("反应时(s)")+ylabel("Density")+ theme_minimal()
