#### STEP 1 ####

##########################
# ACME Synthetic Control 
#   Selection Process
##########################

library(dplyr)
library(tidyr)
library(glmnet)
library(purrr)
library(ggplot2)
library(tibble)

# Load and prepare data
df <- read.csv("ACMEOrdersData.csv") %>%
  mutate(Week_index = row_number())

# Define candidate and full geo sets
candidate_states <- c("Delaware", "Kansas", "Kentucky", "Tennessee",
                      "Missouri", "Montana", "Arkansas", "Mississippi", "New.Mexico")
all_states <- setdiff(names(df), c("Week", "Week_index"))
df$National_Total <- rowSums(df[, all_states])

# Create train/test split
set.seed(42)
train_weeks <- sample(df$Week_index, size = round(0.7 * nrow(df)))
test_weeks <- setdiff(df$Week_index, train_weeks)

# Inputs for LASSO (Step 1)
X_treat <- as.matrix(df[, candidate_states])
y_nat <- df$National_Total

# LASSO model
lasso_model <- glmnet(X_treat[train_weeks, ], y_nat[train_weeks], alpha = 1)

# Define control pool (outside loop)
control_pool <- setdiff(all_states, candidate_states)
control_pool <- control_pool[
  sapply(df[control_pool], function(x) is.numeric(x) && sd(x, na.rm = TRUE) > 1e-5)
]

# Step 1 + Step 2 evaluation over all lambdas
results <- map_dfr(seq_along(lasso_model$lambda), function(i) {
  lambda <- lasso_model$lambda[i]
  coefs <- coef(lasso_model, s = lambda)
  
  selected <- rownames(coefs)[which(coefs != 0)][-1]  # remove intercept
  if (length(selected) < 2 || !all(selected %in% names(df))) return(NULL)
  
  weights <- as.numeric(coefs[which(coefs != 0)][-1])
  weights <- weights / sum(weights)
  
  TG_pred <- as.numeric(as.matrix(df[, selected, drop = FALSE]) %*% weights)
  
  # STEP 1 — How well does TG replicate the National Total?
  step1_mse <- mean((TG_pred[test_weeks] - df$National_Total[test_weeks])^2)
  step1_r2 <- summary(lm(df$National_Total[test_weeks] ~ TG_pred[test_weeks]))$r.squared
  
  # STEP 2 — Use control states to model TG_pred
  control_df <- df %>%
    select(all_of(control_pool)) %>%
    mutate(TG_pred = TG_pred)
  
  model_control <- tryCatch(
    lm(TG_pred ~ ., data = control_df[train_weeks, ]),
    error = function(e) return(NULL)
  )
  if (is.null(model_control)) return(NULL)
  
  synthetic_pred <- predict(model_control, newdata = control_df)
  
  step2_mse <- mean((TG_pred[test_weeks] - synthetic_pred[test_weeks])^2)
  step2_r2 <- summary(lm(TG_pred[test_weeks] ~ synthetic_pred[test_weeks]))$r.squared
  
  tibble(
    Lambda = lambda,
    Selected_TG = paste(selected, collapse = ", "),
    TG_Weights = paste(round(weights, 3), collapse = ", "),
    Step1_MSE = step1_mse,
    Step1_R2 = step1_r2,
    Step2_MSE = step2_mse,
    Step2_R2 = step2_r2
  )
}) %>% arrange(Step2_MSE)

# Show top combinations
print(head(results, 5))
write.csv(results, "two_stage_lasso_synthetic_results.csv", row.names = FALSE)

# Final model visualization using best lambda
best <- results[1, ]
final_selected <- strsplit(best$Selected_TG, ", ")[[1]]
final_weights <- as.numeric(strsplit(best$TG_Weights, ", ")[[1]])
df$treatment_avg <- as.numeric(as.matrix(df[, final_selected, drop = FALSE]) %*% final_weights)

# Rebuild synthetic prediction
final_control_df <- df %>%
  select(all_of(control_pool)) %>%
  mutate(TG_pred = df$treatment_avg)

model_final <- lm(TG_pred ~ ., data = final_control_df[train_weeks, ])
df$synthetic_avg <- predict(model_final, newdata = final_control_df)

# Plot: TG_pred vs Synthetic Control
ggplot(df, aes(x = Week_index)) +
  geom_line(aes(y = treatment_avg, color = "Treatment Avg (TG_pred)"), linewidth = 1.2) +
  geom_line(aes(y = synthetic_avg, color = "Synthetic Control (from Ci)"), linetype = "dashed", linewidth = 1.2) +
  labs(title = paste("Final Synthetic Control vs TG:\n", best$Selected_TG),
       x = "Week", y = "Orders", color = "Legend") +
  scale_color_manual(values = c("Treatment Avg (TG_pred)" = "blue", "Synthetic Control (from Ci)" = "red")) +
  theme_minimal()


### STEP 2 ###
##########################
#    Evaluation of 
#    Campaign Impact 
#   (National Scale)
###########################

# Treatment States (8): Delaware, Kansas, Tennessee, Missouri, Montana, Arkansas, Mississippi, New Mexico

# From previous step:
# final_selected <- vector of selected treatment states
# final_weights <- corresponding weights
# model_final <- pre-trained synthetic control model
# control_pool <- non-candidate states
# df <- historical data frame (Weeks 1–104)


# Load treatment period data
treatment_data <- read.csv("ACME - TreatmentPeriodOrders.csv")
last_week_index <- max(df$Week_index)
treatment_data$Week_index <- seq(last_week_index + 1, last_week_index + nrow(treatment_data))
treatment_data$National_Total <- NA

# Combine with historical data
df_full <- bind_rows(df, treatment_data)

# Define pre- and post-treatment periods
pre_weeks <- df_full$Week_index <= last_week_index
treatment_weeks <- df_full$Week_index > last_week_index

# Stage 1: Predict National_Total from treatment states (holdout – no campaign)
model_stage1 <- lm(National_Total ~ ., data = df_full[pre_weeks, c(final_selected, "National_Total")])
df_full$synthetically_observed_natl <- predict(model_stage1, newdata = df_full)

# Stage 2: Predict National_Total from control states (with campaign)
model_stage2 <- lm(National_Total ~ ., data = df_full[pre_weeks, c(control_pool, "National_Total")])
df_full$synthetic_nat_from_controls <- predict(model_stage2, newdata = df_full)

# Weekly treatment effect: (No Campaign) - (With Campaign)
df_full$national_effect_weekly <- df_full$synthetically_observed_natl - df_full$synthetic_nat_from_controls

# Weekly breakdown
weekly_lift <- df_full[treatment_weeks, c("Week", "Week_index", "synthetically_observed_natl", "synthetic_nat_from_controls", "national_effect_weekly")]
write.csv(weekly_lift, "national_weekly_effects.csv", row.names = FALSE)
print(weekly_lift)

# Total + Average impact
avg_weekly_effect <- mean(df_full$national_effect_weekly[treatment_weeks], na.rm = TRUE)
total_lift <- sum(df_full$national_effect_weekly[treatment_weeks], na.rm = TRUE)
cat("Corrected National Avg Weekly Effect:", round(avg_weekly_effect, 2), "orders/week\n")
cat("Total National Lift (5 Weeks):", round(total_lift, 2), "orders\n")

# Business impact
revenue_lift <- total_lift * 200
profit_lift <- revenue_lift * 0.3
cat("Estimated Revenue Lift: $", round(revenue_lift, 2), "\n")
cat("Estimated Profit Lift: $", round(profit_lift, 2), "\n")

# Load necessary package
library(scales)

# Visualization: Holdout projection vs Synthetic control (National Scale)
ggplot(df_full, aes(x = Week_index)) +
  geom_ribbon(data = df_full[treatment_weeks, ],
              aes(ymin = pmin(synthetically_observed_natl, synthetic_nat_from_controls),
                  ymax = pmax(synthetically_observed_natl, synthetic_nat_from_controls)),
              fill = "lightgray", alpha = 0.4) +
  geom_line(aes(y = synthetically_observed_natl, color = "Observed (No Campaign – Holdout)"), linewidth = 1.2) +
  geom_line(aes(y = synthetic_nat_from_controls, color = "Synthetic Control (With Campaign)"), linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = last_week_index, linetype = "dotted", color = "black") +
  annotate("text", x = last_week_index + 0.5, 
           y = max(df_full$synthetically_observed_natl[treatment_weeks], na.rm = TRUE) * 0.98,
           label = "Campaign Start", angle = 90, vjust = -0.5, size = 3.5) +
  labs(title = "National Causal Impact of Campaign",
       subtitle = "Shaded area = (No Campaign) – (With Campaign)",
       x = "Week Index", y = "National Orders", color = "Legend") +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("Observed (No Campaign – Holdout)" = "blue", 
                                "Synthetic Control (With Campaign)" = "red")) +
  theme_minimal()


# Visualization: Observed vs Counterfactual (Treatment vs Synthetic)
ggplot(df_full[treatment_weeks, ], aes(x = Week_index)) +
  geom_line(aes(y = synthetically_observed_natl, color = "Observed (Treatment Markets)"), linewidth = 1.2) +
  geom_line(aes(y = synthetic_nat_from_controls, color = "Synthetic Control (Counterfactual)"), linetype = "dashed", linewidth = 1.2) +
  labs(title = "Performance Max Campaign Impact (Treatment Period)",
       x = "Week", y = "Orders", color = "Legend") +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("Observed (Treatment Markets)" = "blue", 
                                "Synthetic Control (Counterfactual)" = "red")) +
  theme_minimal()

##########################
#    Null Distribution 
#    via Pre-Treatment 
#    NATIONAL SCALE
##########################

# Null: Weekly differences in pre-period from two national-level synthetic models
df_full$national_null_diff <- df_full$synthetic_nat_from_controls - df_full$synthetically_observed_natl

# Subset to pre-treatment
null_differences <- df_full$national_null_diff[pre_weeks]
observed_effect <- avg_weekly_effect

# Plot aligned to national scale
ggplot(data.frame(diff = null_differences), aes(x = diff)) +
  geom_histogram(fill = "lightgray", color = "black", bins = 30) +
  geom_vline(xintercept = observed_effect, color = "red", linetype = "dashed", linewidth = 1.2) +
  labs(
    title = "Null Distribution of Weekly Differences",
    subtitle = "Pre-Period Weekly Differences vs Observed Campaign Effect",
    x = "Weekly Difference (Observed - Counterfactual)", 
    y = "Frequency"
  ) +
  scale_x_continuous(labels = comma) +
  annotate(
    "text",
    x = observed_effect, 
    y = max(table(cut(null_differences, 30))),
    label = "Observed Effect",
    vjust = -1, hjust = -0.1,
    color = "red"
  ) +
  theme_minimal()


# Compute empirical p-value
p_value <- mean(abs(null_differences) >= abs(observed_effect))
cat("Empirical p-value:", round(p_value, 4), "\n")
