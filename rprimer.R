knitr::opts_chunk$set(echo = TRUE)

## getwd()

## setwd("path/to/directory")

## install.packages("lattice")

## library(lattice)

## # Check current working directory
## getwd()

2 + 3 # Addition
7 - 4 # Subtraction
3 + 5 * 2 # Multiplication
(3 + 5) * 2 # Multiplication with brackets 
#Maths ops follow the normal order of operations.
7/3 # Division
5^2 # Power

150 %/% 60 ## E.g. whole hours in 150 minutes

150 %% 60 ## Minutes left over

x <- 3
x

a <- 1:10 # it increments by one
a

z <- 10

ls()  #list all variables

rm(x)  # delete a variable

ls()

rm(x,y)

rm(list = ls())

## height <-
## 
## weight <-

## bmi <-

#logical
v_log <- c(TRUE, FALSE, FALSE, TRUE) 
v_log

#integer
v_int <- 1:4 
v_int

#numeric double precision
v_doub <- 1:4 * 3.14
v_doub

#character
(v_char <- c("a","b","c","d"))
#note -  I enclosed in brackets to print on screen

is.numeric(v_doub)
is.integer(v_doub)

class(v_doub)

v_log
as.integer(v_log)
v_int
as.numeric(v_int)
v_doub
as.character(v_doub)
as.character(as.numeric(as.integer(v_log)))

x <- c(1.7, "a") ## numeric + character
class(x)
y <- c(TRUE, 2) ## logical + numeric
class(y)
z <- c("a", TRUE) ## character + logical
class(z)

1/0

0/0

v_log <- c(TRUE, FALSE, FALSE, TRUE) #logical
v_int <- 1:4 #integer
v_doub <- 1:4 * 3.14 #numeric double precision
v_char <- c("a","b","c","d") #character

v_char[c(FALSE, FALSE, TRUE, TRUE)]

v_char[v_log]

v_doub[2:3]

v_char[-4]

v_char <- c("Medium","High","Low") # character vector
v_char
v_fac <- factor(v_char) # factor variable
v_fac

levels(v_fac)

class(v_fac)

unclass(v_fac) # In alphabetical order: "Medium=>3,High=>1,Low=>2

v_fac_ord <- factor(v_char, levels = c("Low", "Medium", "High"), ordered = TRUE)
#v_fac_ord <- ordered(v_char, levels = c("Low", "Medium", "High"))
v_fac_ord
unclass(v_fac_ord) 

v_indicator <- c(1, 0, 1, 1, 0, 0,1)

v_fct <- factor(v_indicator,
                levels = c(0, 1),
                labels = c('No', 'Yes'))

summary(v_fct)

f <- factor(c(1, -9, 2))
as.numeric(f)

levels(f)[f] 

as.numeric(levels(f)[f])

as.numeric(as.character(f))

v_fac_ord <- v_fac_ord[v_fac_ord != "Medium"]
summary(v_fac_ord)

f_var_new<-droplevels(v_fac_ord)
summary(f_var_new)

m <- matrix(1:12, nrow = 3, ncol = 4)
m

dim(m)
#Same as:
attributes(m)

my_list <- list(v_char = c('A', 'B'),
                v_num = c(1, 3.5, 8, 15.91),
                v_comp = 1 +4i,
                v_list = list(c(1L, 2L), c(TRUE, FALSE, FALSE)))

my_list

data(iris)

head(iris)

tail(iris)

dim(iris)

nrow(iris)

ncol(iris)

str(iris)

names(iris)

x <- c(1, 2, NA, 10, 3)

is.na(x) # returns logical

is.nan(x)

x <- c(1, 2, NaN, NA, 4)

is.na(x)

is.nan(x)

big <- c(9, 12, 15, 25)
small <- c(9, 3, 4, 2)

big > small

big = small

print(big)
print(small)

big <- c(9, 12, 15, 25)

big[big == small]

big[big > small]

big[big < small]  # Returns an empty set

small>=9 & big<14

big[small>=9 & big<14]

#|
small>=9 | big<14

big[small>=9 | big<14]

x <- c("a", NA, 4, 9, 8.7)
!is.na(x)  # Returns TRUE for non-NA

class(x)

x1 <- x[!is.na(x)]
x1

class(x1)

## data()

## data(package="datasets")

data(iris)

load("./data/demo_data.RData")
str(my.df)

# save specific objects to a file
save(my.df, iris, file = "./data/demo_data1.RData")

# save the whole workspace to the file .RData in the current working directory (cwd)
#save.image()

str(read.csv)

my.df <- read.csv(file = "./data/demo_data.csv")
head(my.df)

## write.csv(my.df, "./data/demo_data_new.csv")

library(readr)
my_df <- read_csv(file = "./data/demo_data.csv")
head(my_df)

## write_csv(my_df, "./data/demo_data_new.csv")

## library(haven)
## #import
## stata_df <- read_dta("./data/demo_data.dta")
## #export
## write_dta(stata_df, "./data/demo_data.dta")

## #import
## sas_df <- read_sas("./data/demo_data.sas7bdat")
## #export
## write_sas(sas_df, "./data/demo_data.sas7bdat")

## #import
## spss_df <- read_sav("./data/demo_data.sav")
## #export
## write_sav(spss_df, "./data/demo_data.sav")

my_df<-read.csv("./data/demo_data.csv")

head(my_df, n=10L)

str(my_df)

summary(my_df)

names(my_df)

#remember $
mean(my_df$age)

mean(my_df$inc1)

mean(my_df$inc1, na.rm=TRUE)

sd(my_df$age) #standard deviation
var(my_df$age) #variance
min(my_df$age) #minimum
max(my_df$age) #maximum
median(my_df$age) #median
range(my_df$age) #range
quantile(my_df$age) #quantile

myvars <- c("inc1", "inc2")
newdata <- my_df[myvars]
#newdata <- my_df[c("inc1","inc2")]
head(newdata,5)

newdata <- my_df[c(3,6:7)]
head(newdata,5)

myvars <- names(my_df) %in% c("inc1", "inc2") 
newdata <- my_df[!myvars]
head(newdata,5)

newdata <- my_df[c(-6,-7)]
head(newdata, n = 5)

#my_df$inc1 <- NULL

newdata <- my_df[1:5,]
newdata

newdata <- my_df[my_df$grp=='C'& my_df$age > 65, ]
head(newdata,5)

str(apply)

apply(my_df[,6:7], 1, sum, na.rm=TRUE)

my_df$tot_inc<-apply(my_df[,6:7], 1, sum, na.rm=TRUE)
head(my_df, 12L)

apply(my_df[,6:7], 2, mean, na.rm=T)

lapply(my_df[6:7], function(x) round(x,1))
#my_df[6:7] <- lapply(my_df[6:7], function(x) round(x,1))

str(tapply)

tapply(my_df$tot_inc, my_df$hhid, sum, na.rm=T)

## install.packages("dplyr")

library("dplyr")

## install.packages("tidyverse")

library("tidyverse")

cat(readLines("./data/nids_subset_description.txt"), sep = '\n')

nids_df <- read_csv("./data/nids_subset.csv")

glimpse(nids_df)

arrange(nids_df, w1_hhid, pid)

filter(nids_df, w1_best_age_yrs>100)

filter(nids_df, w1_best_age_yrs>100 & w1_r_relhead==1)

filter(nids_df, w1_best_age_yrs>100, w1_r_relhead==1)

filter(nids_df, w1_best_age_yrs>100 | w1_r_relhead==1)

select(nids_df,w1_hhid,pid)

nids_df %>%
  filter(w1_r_relhead==1 & w1_best_age_yrs>0) %>%  #filter
  arrange(w1_best_age_yrs) %>% #arrange
  select(w1_hhid, pid, w1_best_gen, w1_best_age_yrs)  #select

## x %>% f %>% g %>% h

## h(g(f(x)))

nids_hh <- nids_df[nids_df$w1_r_relhead==1 & nids_df$w1_best_age_yrs>0, c('w1_hhid', 'pid','w1_best_gen','w1_best_age_yrs')]
nids_hh[order(nids_hh$w1_best_age_yrs),]

nids_df <- nids_df %>% 
  group_by(w1_hhid) %>% 
  arrange(w1_hhid, pid) %>% 
  mutate(pnum = row_number()) #1:n()

nids_df <- nids_df %>% 
  mutate(w1_prov2011 = factor(w1_prov2011, 
                            levels = 1:9,
                            labels = c("Western Cape","Eastern Cape","Northern Cape",
                                       "Free State", "KwaZulu-Natal", "North West",
                                       "Gauteng","Mpumalanga","Limpopo"))) %>% 
  mutate(w1_best_race = factor(w1_best_race,
                               levels = 1:4,
                               c("African","Coloured","Asian_Indian", "White"))) %>% 
  mutate(w1_best_gen = factor(w1_best_gen, levels = 1:2, c("Male","Female")))

glimpse(nids_df)

hhinc_by_prov <- nids_df %>% #create a new data frame from nids_df
  filter(pnum==1) %>% #filter to keep data at the level of the first HH member
  group_by(w1_prov2011) %>% #group by the province variable
  summarize(mean_hhinc = mean(w1_hhincome, na.rm=TRUE)) # compute the mean of hhincome

hhinc_by_prov # print the new data frame

vars_by_prov <- nids_df %>%
  filter(pnum==1) %>% #filter to keep data at the level of the first HH member
    group_by(w1_prov2011) %>% #group by province
    summarize(mean_hhinc = mean(w1_hhincome, na.rm=TRUE), 
              sd_hhinc = sd(w1_hhincome, na.rm=TRUE),
              mean_hhexp = mean(w1_expenditure, na.rm=TRUE),
              sd_hhexp = sd(w1_expenditure, na.rm=TRUE))

vars_by_prov

nids_df %>% 
  filter(pnum==1) %>% #filter to keep data at the level of the first HH member
  group_by(w1_prov2011) %>% #group by province
  summarise(across(c(w1_hhincome,w1_expenditure), list(mean = mean, sd = sd), na.rm=TRUE))

incexppc_byrace_bygen <- nids_df %>%
  group_by(w1_hhid) %>% #group by hhid
  mutate(hhsize = n()) %>% # count number of members
  filter(w1_r_relhead==1) %>% #keep info for HHhead
  mutate(hhinc_pc = w1_hhincome/hhsize, hhexp_pc = w1_expenditure/hhsize) %>% # scale pc
  group_by(w1_best_race, w1_best_gen) %>% #group by race and gender
  summarise(across(c(hhinc_pc,hhexp_pc), list(mean = mean, sd = sd), na.rm=TRUE))

incexppc_byrace_bygen

nids_df %>%
  group_by(w1_prov2011) %>% 
  count(sort = TRUE)

nids_df <- nids_df %>% 
  mutate(wkgrp = if_else(w1_best_age_yrs>=15 & w1_best_age_yrs<=65,1,0),
         dep = if_else(w1_best_age_yrs<15 | w1_best_age_yrs>65,1,0))

glimpse(nids_df)

nids_df<-nids_df %>% 
  mutate(age_bins=case_when(w1_best_age_yrs >= 20 & w1_best_age_yrs <=29 ~ 1,
                             w1_best_age_yrs > 29 & w1_best_age_yrs <=39 ~ 2,
                             w1_best_age_yrs > 39 & w1_best_age_yrs <=49 ~ 3,
                             w1_best_age_yrs > 49 & w1_best_age_yrs <=59 ~ 4,
                             w1_best_age_yrs > 59 & w1_best_age_yrs <=69 ~ 5,
                             w1_best_age_yrs > 69 & w1_best_age_yrs <=120 ~ 6))

nids_df <- nids_df %>% 
  mutate(age_bins = factor(age_bins,
                           levels = 1:6,
                           labels = c("20 - 29 yrs","30 - 39 yrs", "40 - 49 yrs", "50 - 59 yrs",
                                      "60 - 69 yrs", "70 - 120 yrs")))

nids_df%>%
  group_by(age_bins)%>%
  summarise(freq = n())

#nids_df %>%  group_by(age_bins) %>% count

inc_exp <- nids_df %>% 
  filter(w1_r_relhead==1) %>% 
  select(w1_hhid, w1_hhincome,w1_expenditure, w1_best_race, w1_best_gen)

str(plot)

plot(inc_exp$w1_hhincome)

plot(inc_exp$w1_hhincome, inc_exp$w1_expenditure)

plot(log10(inc_exp$w1_hhincome),log10(inc_exp$w1_expenditure), 
     xlab = "log(Household income)", #horizontal axis label
     ylab = "log(Household expenditure)", #vertical axis label
     main = "Variation of household income and expenditure", #plot title
     col = "blue", #colour
     pch = 20) #sympol

hist(inc_exp$w1_hhincome)

hist(inc_exp$w1_hhincome[inc_exp$w1_hhincome<=50000], 
     breaks=50,
     xlab = "Household income",
     main = "Distribution of household income from NIDS subset")

counts<-table(nids_df$w1_best_race)
barplot(counts)

counts1<-table(nids_df$w1_best_gen, nids_df$w1_best_race)

barplot(counts1, 
        main="Population group by gender",
        xlab="Population group", 
        #col=c("blue","red"),
        legend = rownames(counts1))

barplot(counts1, 
        main="Population group by gender",
        xlab="Population group", 
        col=c("blue","red"),
        legend = rownames(counts1), beside=TRUE)

# we excluded expenditure greater than 15k for ease of viewing the plot 
boxplot(inc_exp$w1_expenditure[inc_exp$w1_expenditure<15000]) 

boxplot(w1_expenditure ~ w1_best_race, 
        data=inc_exp[inc_exp$w1_expenditure<15000,],
        main = "Boxplot of household expenditure",
        xlab = "Population group", ylab = "Household expenditure",cex.lab = 1.25)

plot(inc_exp$w1_hhincome,inc_exp$w1_expenditure,
     xlab = "Household income", #horizontal axis label
     ylab = "Household expenditure", #vertical axis label
     main = "Variation of household income and expenditure", #plot title
     col = "blue", #colour
     pch = 20) #sympol
abline(a = 1137, b = 0.61,lwd = 2, col="red") # add abline

plot(inc_exp$w1_hhincome,inc_exp$w1_expenditure, 
     xlab = "Household income", #horizontal axis label
     ylab = "Household expenditure", #vertical axis label
     main = "Variation of household income and expenditure", #plot title
     col = "blue", #colour
     pch = 20) #sympol
points(inc_exp[inc_exp$w1_expenditure>inc_exp$w1_hhincome, c("w1_hhincome", "w1_expenditure")], col = "red")

plot(inc_exp[c("w1_hhincome","w1_expenditure")], 
     xlab = "Household income", #horizontal axis label
     ylab = "Household expenditure", #vertical axis label
     main = "Variation of household income and expenditure", #plot title
     pch = 20) #sympol
lines(stats::lowess(inc_exp[c("w1_hhincome","w1_expenditure")]), col = "blue")

plot(inc_exp$w1_hhincome,inc_exp$w1_expenditure, 
     xlab = "HH income", #horizontal axis label
     ylab = "HH expenditure", #vertical axis label
     main = "Variation of household income and expenditure", #plot title
     col = "blue", #colour
     pch = 20) #sympol
abline(a = 1137, b = 0.61, lwd = 2, col = "red") #add line
text(x = 100000,y = 75000, labels = "y = 1137 + 0.61x") #add equation at (x,y)

pdf("./figures/myplot.pdf")
boxplot(nids_df$w1_hhincome ~ nids_df$w1_best_race, main = "Boxplot of HH income",
        xlab = "Population group", ylab = "Household income",cex.lab = 1.25)
dev.off()

png(filename = "./figures/myplot.png", width = 20, height = 20, units = "cm", res = 400)
boxplot(nids_df$w1_hhincome ~ nids_df$w1_best_race, main = "Boxplot of HH income",
        xlab = "Population group", ylab = "Household income",cex.lab = 1.25)
dev.off()

library(ggplot2)

ggplot(data = inc_exp)

ggplot(data = inc_exp, aes(x=w1_hhincome, y=w1_expenditure))

ggplot(data = inc_exp, aes(x=w1_hhincome, y=w1_expenditure)) + 
  geom_point()

p <- ggplot(data = inc_exp, aes(x = w1_hhincome, y = w1_expenditure)) +
  geom_point(alpha=0.3) +  # transparency
  geom_smooth(method="lm") + # smooth distribution - linear model
  scale_x_log10(labels = scales::dollar_format(suffix = "", prefix = "R")) + #log scale
  scale_y_log10(labels = scales::dollar_format(suffix = "", prefix = "R")) + #log scale
  labs(title = "Household expenditure by income", #title
       subtitle = "From NIDS 2008", # subtitle
       x = "Household income", # x-axis title
       y = "Household expenditure", # y-axis title
       caption="Source: NIDS") + # caption
  theme_minimal() # change theme
p

ggplot(inc_exp, aes(w1_hhincome))+
  geom_histogram()

ggplot(inc_exp %>% filter(w1_hhincome<50000), aes(w1_hhincome))+
  geom_histogram(binwidth = 750)

ggplot(inc_exp %>% filter(w1_hhincome<50000), aes(w1_hhincome))+
  geom_histogram()+
  scale_x_log10()

inc_exp %>% filter(w1_expenditure<15000) %>% 
  ggplot(aes(x = w1_best_race, y = w1_expenditure))+
  geom_boxplot()

ggplot(data = inc_exp, aes(x=w1_hhincome, y=w1_expenditure, shape = w1_best_race)) + 
  geom_point() 

ggplot(data = inc_exp, aes(x=w1_hhincome, y=w1_expenditure, color = w1_best_race )) + 
  geom_point() 

ggplot(inc_exp, aes(x = w1_hhincome, y=w1_expenditure)) + 
  geom_point() + 
  geom_smooth() +
  facet_wrap(~w1_best_race)

#save last displayed plot as pdf
ggsave(plot = p, "./figures/incVexp.pdf")
#save plot stored in "p" as png and resize
ggsave(plot = p, "./figures/incVexp.png", width=10, height=10, units="in")
#?ggsave

xx <- data.frame(x = c(1,2.5,3), y = c(2,3,1))

theme_noaxis <- theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank())

p1 <- ggplot(xx, aes(x,y))+
  geom_point()+
  labs(x = "", y = "", title = "Points")+
  coord_fixed(ratio = 1)+
  theme_minimal()+
  theme_noaxis

p2 <- ggplot(xx, aes(x,y))+
  geom_point()+
  geom_line()+
  labs(x = "", y = "", title = "Line")+
  coord_fixed(ratio = 1)+
  theme_minimal()+
  theme_noaxis

p3 <- ggplot(xx, aes(x,y))+
  geom_point()+
  geom_polygon(alpha=0.3)+
  labs(x = "", y = "", title = "Polygon")+
  coord_fixed(ratio = 1)+
  theme_minimal()+
  theme_noaxis

library(patchwork)
p1 + p2 + p3

# load packages
library(sf)

sa_prov<-st_read(dsn="./PR_SA_2011", layer="PR_SA_2011") #no extension
#sa_prov<-st_read(dsn="./shapefiles/PR_SA_2011/PR_SA_2011.shp") #extension

class(sa_prov)

sa_prov

# pull just the geometry
sa_prov %>% st_geometry() %>% plot()
#plot(st_geometry(sa_prov))

plot(sa_prov["SHAPE_Leng"])

sessionInfo()
