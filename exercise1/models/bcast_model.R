setwd(getwd())

getwd()

data_bcast = read.csv("../results_def/bcast_complete.csv")

# Correlation map
cor(data_bcast)

# Convert to factors the columns Allocation and Algorithm
data_bcast$Allocation = as.factor(data_bcast$Allocation)
data_bcast$Algorithm = as.factor(data_bcast$Algorithm)

head(data_bcast)

# Create a linear model
model = lm(Avg.Latency.us. ~ . , data = data_bcast)

summary(model)

# Merge the default factor levels with bin_tree
# library(forcats)
# data_bcast$Algorithm = fct_collapse(data_bcast$Algorithm, "bin_tree" = c("bin_tree", "default", "linear"))

# Create a linear model
model1 = lm(Avg.Latency.us. ~ Algorithm + Processes + MessageSize , data = data_bcast[data_bcast$Allocation == "core",])

summary(model1)

# Try to take into consideration interaction
model2 = lm(Avg.Latency.us. ~ Algorithm * Allocation + Allocation * Processes + Allocation * Processes * Algorithm , data = data_bcast)

summary(model2)



#try with a GAM
library(mgcv)
model_gam = gam(Avg.Latency.us. ~ Algorithm + Allocation + s(Processes) + MessageSize , data = data_bcast)

summary(model_gam)
