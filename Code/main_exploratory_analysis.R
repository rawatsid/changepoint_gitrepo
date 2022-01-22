library(ggplot2)
library(ggpubr)
library(usmap)
library(spacetime)

###Exploratory analysis#######################################################################################

data_weekly_df = read.csv("D:/Year_3/Term2/Research/changepoint_analysis/US_COVID_weekly_data.csv",row.names = NULL)

data_weekly_df$week = as.Date(data_weekly_df$week)

data_weekly_NY_df = data_weekly_df %>% filter(state == "New York")

data_weekly_NY_df %>% filter(id=="ID1946")  %>% #Richmond
  ggplot(aes(x = week,y = new_cases)) +
  geom_line() +
  #geom_point(aes(x = week,y = new_cases),size = 3) +
  #scale_shape_manual(values = c("significant" = 19,"not significant" = 1)) +
  ylab("Richmond New Cases") +
  xlab("Week")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.title = element_blank())

data_weekly_NY_df %>% filter(id =="ID1946")  %>% #Richmond
  ggplot(aes(x = week,y = new_cases_per_100k)) +
  geom_line() +
  #geom_point(aes(x = week,y = new_cases),size = 3) +
  #scale_shape_manual(values = c("significant" = 19,"not significant" = 1)) +
  ylab("Richmond New Cases per 100k") +
  xlab("Week")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.title = element_blank())

data_weekly_NY_df %>% filter(id=="ID1926")  %>% #Kings
  ggplot(aes(x = week,y = new_cases)) +
  geom_line() +
  #geom_point(aes(x = week,y = new_cases),size = 3) +
  #scale_shape_manual(values = c("significant" = 19,"not significant" = 1)) +
  ylab("Kings New Cases") +
  xlab("Week")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.title = element_blank())

#New state level data
data_weekly_NY_df %>% group_by(country,state,week) %>%summarise(new_cases_state=sum(new_cases)) %>% #Kings
  ggplot(aes(x = week,y = new_cases_state)) +
  geom_line() +
  #geom_point(aes(x = week,y = new_cases),size = 3) +
  #scale_shape_manual(values = c("significant" = 19,"not significant" = 1)) +
  ylab("New York New Cases") +
  xlab("Week")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.title = element_blank())

fips_code_df = read.csv("D:/Year_3/Term2/Research/changepoint_analysis//FIPS_code_NY.csv", row.names = NULL)

data_weekly_NY_df = data_weekly_NY_df %>% 
  mutate(county = gsub("^(.*?),.*", "\\1", key))

names(fips_code_df)[1] = "fips"
names(fips_code_df)[2] = "county"

fips_code_df[45,2]= "St. Lawrence"

data_weekly_NY_fips_df = merge(fips_code_df[,1:2],data_weekly_NY_df, by = "county")

data_weekly_NY_fips_df = data_weekly_NY_fips_df %>% 
  arrange(country, state, county,lat, long, id, key, week)



#Get fips code

# divide the data according to quartiles
data_subset1 = data_weekly_NY_fips_df  %>%
  filter(week == as.Date("2020-06-29")) %>%
  dplyr::select(fips,confirmed) 

data_subset2 = data_weekly_NY_fips_df  %>%
  filter(week == as.Date("2020-12-21")) %>%
  dplyr::select(fips,confirmed) 

data_subset3 = data_weekly_NY_fips_df  %>%
  filter(week == as.Date("2021-06-07")) %>%
  dplyr::select(fips,confirmed) 

data_subset4 = data_weekly_NY_fips_df  %>%
  filter(week == as.Date("2021-11-22")) %>%
  dplyr::select(fips,confirmed) 

# set lower and upper limit of every plot
LL = min(data_weekly_NY_df$new_cases)
UL = max(data_weekly_NY_df$new_cases)

# create four different plots
p1 = plot_usmap(data = data_subset1,values = "confirmed",color = "black",include = c("NY")) + 
  scale_fill_continuous(name = "confirmed",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("29 Jun, 2020") +
  theme(plot.title = element_text(hjust = 0.5,size = 35),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))

p2 = plot_usmap(data = data_subset2,values = "confirmed",color = "black",include = c("NY")) + 
  scale_fill_continuous(name = "confirmed",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("21 Dec, 2020") +
  theme(plot.title = element_text(hjust = 0.5,size = 35),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))

p3 = plot_usmap(data = data_subset3,values = "confirmed",color = "black",include = c("NY")) + 
  scale_fill_continuous(name = "confirmed",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("7 Jun, 2021") +
  theme(plot.title = element_text(hjust = 0.5,size = 35),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))

p4 = plot_usmap(data = data_subset4,values = "confirmed",color = "black",include = c("NY")) + 
  scale_fill_continuous(name = "confirmed",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("22 Nov, 2021") +
  theme(plot.title = element_text(hjust = 0.5,size = 35),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))

# merge them together and save on file
ggpubr::ggarrange(p1,p2,p3,p4)
######################################################################################################
#Proportion of time category
data_weekly_NY_fips_df$lat = round(data_weekly_NY_fips_df$lat, digits = 3)
data_weekly_NY_fips_df$long = round(data_weekly_NY_fips_df$long, digits = 3)


data_weekly_NY_fips_countycateg_df = data_weekly_NY_fips_df %>% 
  group_by(country, state, county,lat, long, id, key, fips, population) %>% summarise(
    prop_categ1 = sum(category == 1)/n(),
    prop_categ2 = sum(category == 2)/n(),
    prop_categ3 = sum(category == 3)/n(),
    prop_categ4 = sum(category == 4)/n())

data_subset1 = data_weekly_NY_fips_countycateg_df  %>%
  dplyr::select(fips,prop_categ1) 

data_subset2 = data_weekly_NY_fips_countycateg_df  %>%
  dplyr::select(fips,prop_categ2) 

data_subset3 = data_weekly_NY_fips_countycateg_df  %>%
  dplyr::select(fips,prop_categ3) 

data_subset4 = data_weekly_NY_fips_countycateg_df  %>%
  dplyr::select(fips,prop_categ4) 

# set lower and upper limit of every plot
LL = min(data_weekly_NY_fips_countycateg_df$prop_categ1,data_weekly_NY_fips_countycateg_df$prop_categ2,
         data_weekly_NY_fips_countycateg_df$prop_categ3, data_weekly_NY_fips_countycateg_df$prop_categ4)
UL = max(data_weekly_NY_fips_countycateg_df$prop_categ1,data_weekly_NY_fips_countycateg_df$prop_categ2,
         data_weekly_NY_fips_countycateg_df$prop_categ3, data_weekly_NY_fips_countycateg_df$prop_categ4)

# create four different plots
p1 = plot_usmap(data = data_subset1,values = "prop_categ1",color = "black",include = c("NY")) + 
  scale_fill_continuous(name = "prop_categ1",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("Proportional Category 1") +
  theme(plot.title = element_text(hjust = 0.5,size = 35),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))

p2 = plot_usmap(data = data_subset2,values = "prop_categ2",color = "black",include = c("NY")) + 
  scale_fill_continuous(name = "prop_categ2",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("Proportional Category 2") +
  theme(plot.title = element_text(hjust = 0.5,size = 35),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))

p3 = plot_usmap(data = data_subset3,values = "prop_categ3",color = "black",include = c("NY")) + 
  scale_fill_continuous(name = "prop_categ3",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("Proportional Category 3") +
  theme(plot.title = element_text(hjust = 0.5,size = 35),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))

p4 = plot_usmap(data = data_subset4,values = "prop_categ4",color = "black",include = c("NY")) + 
  scale_fill_continuous(name = "prop_categ4",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("Proportional Category 4") +
  theme(plot.title = element_text(hjust = 0.5,size = 35),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))

# merge them together and save on file
ggpubr::ggarrange(p1,p2,p3,p4)
#######################################################################################################
#Location-wise increasing or decreasing direction##############
data_weekly_NY_fips_categ_inc_dec_df = data_weekly_NY_fips_df %>% 
  group_by(country, state, county,lat, long, id, key, fips, population) %>% summarise(
    count_inc_categ = sum(diff(category) > 0),
    count_dec_categ = sum(diff(category) < 0)
  )

data_subset1 = data_weekly_NY_fips_categ_inc_dec_df  %>%
  dplyr::select(fips,count_inc_categ) 

data_subset2 = data_weekly_NY_fips_categ_inc_dec_df  %>%
  dplyr::select(fips,count_dec_categ) 

 

# set lower and upper limit of every plot
LL = min(data_weekly_NY_fips_categ_inc_dec_df$count_inc_categ, 
         data_weekly_NY_fips_categ_inc_dec_df$count_dec_categ)
UL = max(data_weekly_NY_fips_categ_inc_dec_df$count_inc_categ,
         data_weekly_NY_fips_categ_inc_dec_df$count_dec_categ)

# create four different plots
p1 = plot_usmap(data = data_subset1,values = "count_inc_categ",color = "black",include = c("NY")) + 
  scale_fill_continuous(name = "count_inc_categ",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("Increasing Number of Categories") +
  theme(plot.title = element_text(hjust = 0.5,size = 35),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))

p2 = plot_usmap(data = data_subset2,values = "count_dec_categ",color = "black",include = c("NY")) + 
  scale_fill_continuous(name = "count_dec_categ",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("Decreasing Number of Categories") +
  theme(plot.title = element_text(hjust = 0.5,size = 35),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))



# merge them together and save on file
ggpubr::ggarrange(p1,p2,ncol = 1)
#######################################################################################################
#Spatial autocorrelation#####################################
library(spdep)
library(sf)

spdf <- training_df
coordinates(spdf) <- ~ long + lat

# create the neighbours list for the SLX code
neib <- knn2nb(knearneigh(coordinates(spdf),longlat = TRUE))
lw <- nb2listw(neib,style = "B")

joincount.test(as.factor(training_df$category),lw)



