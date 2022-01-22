library(tidyverse)
library(dplyr)


data_df = read.csv("D:/Year_3/Term2/Research/changepoint_analysis/US_COVID_data.csv", row.names = NULL)
data_df$date = as.Date(data_df$date)
# add week column (last Monday) to the data
lastmon2 <- function(x) x - as.numeric(x-1+4)%%7

data_df$lat = round(data_df$lat, digits = 3)
data_df$long = round(data_df$long, digits = 3)

data_df = data_df %>%
  mutate(
    week = lastmon2(date)
  )
#######To check if id and key are unique one to one mapping###########################
View(unique(data_df %>% select(id, key)) %>% group_by(id)  %>% summarise(count_id = n()))
View(unique(data_df %>% select(id, key)) %>% group_by(key)  %>% summarise(count_id = n()))

data_weekly_df = data_df %>%
  group_by(country,state,lat,long,id,key,population,week) %>%
  summarize(
    confirmed = max(confirmed), # getting total cases as the maximum in that week
    death = max(death)
  ) %>%
  ungroup() %>% arrange(country, state, key, lat, long, id,population,  week) %>% 
  group_by(country,state,key,lat,long,id,population) %>%
  mutate(
    time = c(1:n())/n(),
    time_sq = time^2,
    prev_week_total = dplyr::lag(confirmed,default = min(confirmed)),
    new_cases = confirmed - prev_week_total,
    incidence = ifelse(new_cases > 0,new_cases/population,0),
    log_incidence = log(incidence + 0.1), # avoiding error for 0 cases
    prev_week_actual = dplyr::lag(new_cases,default = min(confirmed)),
    
    log_death = log(death + 0.1),
    prev_log_death = lag(log_death,1,default = log_death[1]),
    new_death = death - lag(death,1,default = 0),
    new_death = ifelse(new_death >= 0,new_death,0),    # correcting one row where new_death<0
    log_new_death = log(new_death + 0.1),
    prev_log_new_death = lag(log_new_death,1,default = log_new_death[1]),
    new_cases_per_100k = new_cases * 100000 / population,
    category = ifelse(new_cases_per_100k>=100, 4, ifelse(new_cases_per_100k>=50, 3, 
                                                         ifelse(new_cases_per_100k>=10,2,1)))
  ) %>%
  ungroup() %>% arrange(country, state, key, lat, long, id,population,  week)


write.csv(x = data_weekly_df,file = "D:/Year_3/Term2/Research/changepoint_analysis/US_COVID_weekly_data.csv",
          row.names = F)

data_weekly_NY_df = data_weekly_df %>% filter(state == "New York")

write.csv(x = data_weekly_NY_df,file = "D:/Year_3/Term2/Research/changepoint_analysis/US_COVID_weekly_NY_data.csv",
          row.names = F)

