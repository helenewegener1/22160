---
title: "Helene Wegener (s165927)"
output: html_document
date: "2023-07-27"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Load libraries  
```{r}
library(tidyverse)
```

# Load data
In terminal: git clone https://github.com/ramhiser/datamicroarray.git
Open "gravier" in environment. 

# Analysis 
## 1. Create a tibble from the gravier data
```{r}
data <- gravier$x %>% as_tibble()
```

## 2. Relocate the “y” column to the first column and rename it “outcome”
```{r}
data <- cbind(outcome = gravier$y, data)
```

## 3. Recode outcome, such that 0 is “good” and 1 is “poor”
```{r}
data <- data %>% 
  mutate(outcome = case_when(outcome == 'good' ~ 0,
                             outcome == 'poor' ~ 1))
```

## 4. That is your clean data, from that...

## 5. Transform to long format
```{r}
data <- data %>% 
  pivot_longer(cols = starts_with('g'),
               names_to = 'gene',
               values_to = 'log2_expr_lvl')
```

## 6. Choose e.g. 100 random genes
```{r}
random_genes <- unique(data$gene) %>% 
  sample(100)
```

## 7. Fit a logistic regression to each gene, modelling: “outcome ~ log2_expr_lvl”
```{r}
logistic_regression <- function(gene){
  
  model <- glm(formula = outcome ~ log2_expr_lvl, family = binomial(link = "logit"), data = data[data$gene == gene, ])
  
  beta <- model$coefficients['log2_expr_lvl']
  names(beta) <- NULL
  
  conf <- confint(model)['log2_expr_lvl', ]
  names(conf) <- NULL
  conf_2.5 <- conf[1]
  conf_97.5 <- conf[2]
  
  pvalue <- coef(summary(model))['log2_expr_lvl', 'Pr(>|z|)']
  
  return(c(beta, conf_2.5, conf_97.5, pvalue))
  
}

fit <- data.frame('gene' = random_genes,
                  'beta-estimates' = sapply(random_genes, logistic_regression)[1 ,],
                  'conf_2.5' = sapply(random_genes, logistic_regression)[2 ,], 
                  'conf_97.5' = sapply(random_genes, logistic_regression)[3 ,],
                  'pvalue' = sapply(random_genes, logistic_regression)[4 ,], 
                  row.names = NULL) 
```

## 8. Add beta-estimates and confidence intervals
Done above.

## 9. Add an indicator for whether the p-value was less than or equal to 0.05
```{r}
fit <- fit %>% 
  mutate(is_significant = case_when(pvalue <= 0.05 ~ 'significant',
                                    pvalue > 0.05 ~ 'n.s'))
```

## 10. That is your long modelled data, from that:

## 11. Create a forest-plot of the slopes (beta1 estimates) and add 95% CI
```{r}
fit %>% 
  mutate(gene = fct_reorder(gene, desc(beta.estimates))) %>% 
  ggplot(aes(x = beta.estimates,
             y = reorder(gene, beta.estimates, decreasing = TRUE),
             color = is_significant)) +
  geom_point() +
  geom_errorbar(aes(xmin = conf_2.5, 
                     xmax = conf_97.5)) +
  geom_vline(xintercept = 0) + 
  labs(title = 'Gravier: Gene expression estiamtes',
       x = 'beta estimates',
       y = 'gene') + 
  theme_light()

ggsave('forestplot.pdf',
       width = 10,
       height = 11)
```
