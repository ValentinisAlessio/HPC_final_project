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

compl_model <-summary(model)
compl_model_table <- xtable(compl_model)
print.xtable(compl_model_table, file = "compl_model_table_bar.tex", floating = FALSE, type = "latex")

# Consider only core allocation
model1 = lm(Avg.Latency.us. ~ Algorithm + Processes , data = data_barrier[data_barrier$Allocation == "core",])

summary(model1)
fixed_model <-summary(model1)
fixed_model_table <- xtable(fixed_model)
print.xtable(fixed_model_table, file = "fixed_model_table_bar.tex", floating = FALSE, type = "latex")

#try with a GAM
library(mgcv)
model_gam_barrier = gam(Avg.Latency.us. ~ Algorithm + Allocation + s(Processes), data = data_barrier)

summary(model_gam_barrier)

