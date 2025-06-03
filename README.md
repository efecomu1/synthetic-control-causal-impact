# Synthetic Control for Geo-Test Campaign Analysis (ACME)

This project evaluates the causal impact of a Google PMAX campaign on retail order volume using Lasso regression and synthetic control methods.

## 📌 Objective
- Estimate national campaign lift using a geo-holdout and synthetic control strategy.
- Quantify causal effect on order volume, revenue, and profit.

## 🧠 Methodology
- **Step 1**: Use Lasso to identify treatment group (TG) from candidate states.
- **Step 2**: Create a synthetic control using other states as predictors.
- **Step 3**: Project counterfactual national trends and compute treatment effect.

## 📊 Key Results
- National weekly lift: +26,448 orders
- Total revenue lift (5 weeks): $5.3M
- Empirical p-value: 0 (strong evidence of causal effect)

## 🛠️ Tools & Tech
- R (glmnet, ggplot2, dplyr)
- Synthetic control, Lasso regression
- Pre/post analysis and null distribution validation

## 📄 Files
- `scripts/causal_impact_analysis.R`: Full analysis code
- `presentation/`: Slide deck from final class presentation
