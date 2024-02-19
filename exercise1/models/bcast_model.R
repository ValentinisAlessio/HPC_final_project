setwd(getwd())

getwd()

library(xtable)

data_bcast = read.csv("../results_def/bcast_complete.csv")

# Correlation map
cor(data_bcast)

# Convert to factors the columns Allocation and Algorithm
data_bcast$Allocation = as.factor(data_bcast$Allocation)
data_bcast$Algorithm = as.factor(data_bcast$Algorithm)

head(data_bcast)

# Create a linear model
model = lm(Avg.Latency.us. ~ . , data = data_bcast)

compl_model <-summary(model)
compl_model_table <- xtable(compl_model)
print.xtable(compl_model_table, file = "compl_model_table.tex", floating = FALSE, type = "latex")


# Merge the default factor levels with bin_tree
# library(forcats)
# data_bcast$Algorithm = fct_collapse(data_bcast$Algorithm, "bin_tree" = c("bin_tree", "default", "linear"))

# Create a linear model
model1 = lm(Avg.Latency.us. ~ Algorithm + Processes + MessageSize , data = data_bcast[data_bcast$Allocation == "core",])

summary(model1)
fixed_model <-summary(model1)
fixed_model_table <- xtable(fixed_model)
print.xtable(fixed_model_table, file = "fixed_model_table.tex", floating = FALSE, type = "latex")

# Try to take into consideration interaction
model2 = lm(Avg.Latency.us. ~ Algorithm + Processes + poly(MessageSize, degree=), data = data_bcast[data_bcast$Allocation == "core",])

summary(model2)



#try with a GAM
library(mgcv)
model_gam = gam(Avg.Latency.us. ~ Algorithm + s(Processes) + s(MessageSize) , data = data_bcast)

compl_gam_model <- summary(model_gam)
compl_gam_model_table <- xtable(compl_gam_model)
print.xtable(compl_gam_model_table, file = "compl_gam_model_table.tex", floating = FALSE, type = "latex")

model_gam = gam(Avg.Latency.us. ~ Algorithm + s(Processes) + s(MessageSize) , data = data_bcast[data_bcast$Allocation == "core",])