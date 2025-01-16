# Install and load the xgboost package
# install.packages("xgboost")
# install.packages("caret")

# Other useful libraries
library(xgboost)  # Machine Learning
library(dplyr)    # Data manipulation
library(caret)    # For creating a train-test split
library(pROC)
library(ggplot2)

# set working directory
setwd("C:/Users/qdani/OneDrive/Documents/grad_school/coursework/ECON8210/replication")

# Set seed for reproducibility
set.seed(123)

# Generate Panel Data

# # Set number of agents
# n = 123000
# 
# # Set default probability of 10% among approved users, independent of gender
# y <- sample(c(rep(0,9),1),n,replace = TRUE)
# 
# # Set gender shares , 1 = male, 0 = female
# gender <- sample(c(0,1),n,replace = TRUE, prob = c(.38,.62))
# 
# #Assemble into panel
# panel <- data.frame(y, # Default outcome
#                     gender,
#                     age=vector(length=n), #Starting here, empty vectors
#                     ios=vector(length=n),
#                     nohit=vector(length=n),
#                     num_ord=vector(length=n),
#                     cash_p=vector(length=n),
#                     med_amt=vector(length=n),
#                     supm_p=vector(length=n),
#                     phrm_p=vector(length=n),
#                     rest_p=vector(length=n),
#                     ses=vector(length=n),
#                     years_ed=vector(length=n),
#                     mveh_p=vector(length=n)
#                     )
# 
# # Approximation -- eyeballed it
# for (i in 1:n) {
#   
#   slice <- panel[i,]
#   default <- slice$y
#   gender <- slice$gender
#   
#              if(default == 1 & gender == 1){
#     slice$ios = sample(c(0,1),size = 1,prob = c(.83,.17))
#     slice$age = max(rnorm(1,23.3,7.3),18)
#     slice$nohit = rnorm(1,630,18)
#     slice$num_ord = rnorm(1,10,3) + runif(1,0,45)
#     slice$cash_p = max(min(rnorm(1,.81,.36),1),0)
#     slice$med_amt = rnorm(1,200,50) + runif(1,0,350)
#     slice$supm_p = rnorm(1,.05,.01) + runif(1,0,.1)
#     slice$phrm_p = rnorm(1,.05,.01) + runif(1,0,.08)
#     slice$rest_p = min(rnorm(1,.83,.27),1)
#     slice$ses = rnorm(1,.96,.01)
#     slice$years_ed = rnorm(1,11,1.7)
#     slice$mveh_p = rnorm(1,.52,.17)
#       } else if(default == 0 & gender == 1){
#     slice$ios = sample(c(0,1),size = 1,prob = c(.64,.36))
#     slice$age = max(rnorm(1,24.8,8.3),18)
#     slice$nohit = rnorm(1,641,20.7)
#     slice$num_ord = rnorm(1,10,3) + runif(1,0,45)
#     slice$cash_p = min(max(rnorm(1,.48,.36),0),1)
#     slice$med_amt = rnorm(1,150,50) + runif(1,0,300)
#     slice$supm_p = rnorm(1,.06,.01) + runif(1,0,.1)
#     slice$phrm_p = rnorm(1,.06,.01) + runif(1,0,.08)
#     slice$rest_p = min(rnorm(1,.79,.27),1)
#     slice$ses = rnorm(1,.97,.01)
#     slice$years_ed = rnorm(1,12.5,1.7)
#     slice$mveh_p = rnorm(1,.64,.17)
#       } else if (default == 1 & gender == 0){
#     slice$ios = sample(c(0,1),size = 1,prob = c(.79,.21))
#     slice$age = max(rnorm(1,24.8,7.3),18)
#     slice$nohit = rnorm(1,630,18)
#     slice$num_ord = rnorm(1,11,3) + runif(1,0,45)
#     slice$cash_p = max(min(rnorm(1,.81,.36),1),0)
#     slice$med_amt = rnorm(1,200,50) + runif(1,0,350)
#     slice$supm_p = rnorm(1,.05,.01) + runif(1,0,.1)
#     slice$phrm_p = rnorm(1,.05,.01) + runif(1,0,.08)
#     slice$rest_p = min(rnorm(1,.83,.27),1)
#     slice$ses = rnorm(1,.96,.01)
#     slice$years_ed = rnorm(1,11,1.7)
#     slice$mveh_p = rnorm(1,.52,.17)
#       } else if (default == 0 & gender == 0){
#     slice$ios = sample(c(0,1),size = 1,prob = c(.58,.42))
#     slice$age = max(rnorm(1,25.8,8.3),18)
#     slice$nohit = rnorm(1,641,20.7)
#     slice$num_ord = rnorm(1,12,3) + runif(1,0,45)
#     slice$cash_p = min(max(rnorm(1,.48,.36),0),1)
#     slice$med_amt = rnorm(1,150,50) + runif(1,0,300)
#     slice$supm_p = rnorm(1,.06,.01) + runif(1,0,.1)
#     slice$phrm_p = rnorm(1,.06,.01) + runif(1,0,.08)
#     slice$rest_p = min(rnorm(1,.79,.27),1)
#     slice$ses = rnorm(1,.97,.01)
#     slice$years_ed = rnorm(1,12.5,1.7)
#     slice$mveh_p = rnorm(1,.64,.17)
#       }
#     
#   
#   panel[i,] <- slice
#   
#   print(i)
#   
# }
# 
# saveRDS(panel,"grad_school/coursework/ECON8210/replication/panel.RDS")
