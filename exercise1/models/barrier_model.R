setwd(getwd())

getwd()

data_barrier = read.csv("../results_def/barrier_complete.csv")

# Convert to factors the columns Allocation and Algorithm
data_barrier$Allocation = as.factor(data_barrier$Allocation)
data_barrier$Algorithm = as.factor(data_barrier$Algorithm)

cor(data_barrier)

head(data_barrier)

# Create a linear model
model_barrier = lm(Avg.Latency.us. ~ Algorithm + Allocation + Processes , data = data_barrier)

summary(model_barrier)

#try with a GAM
library(mgcv)
model_gam_barrier = gam(Avg.Latency.us. ~ Algorithm + Allocation + s(Processes), data = data_barrier)

summary(model_gam_barrier)

