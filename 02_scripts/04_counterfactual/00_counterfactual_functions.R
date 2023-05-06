# Functions for 01_counterfactual_model.Rmd
## Functions for Counterfactual Analysis

### Function for counterfactual model
model_counter = function(data, 
                         u_d = u_d, 
                         u_j = u_j, 
                         iter = 2000, 
                         chains = 4, 
                         seed = 34568979){
  
  # Define protected variables
  data_prot = data %>% 
    select(starts_with(c("Race_", "Sex_"))) 
  
  n = nrow(data)
  
  # Create list for Stan file
  data_list= list(N = n,
                  K = length(data_prot),
                  prot = data.matrix(data_prot),
                  age = data[, "Age", drop = TRUE],
                  priors = data[, "Priors", drop = TRUE],
                  crime = data[, "Crime", drop = TRUE],
                  juvmisd = data[, "JuvMisd", drop = TRUE],
                  juvfel = data[, "JuvFel", drop = TRUE],
                  compas = data[, "COMPAS", drop = TRUE],
                  mu_d = u_d[, "median", drop = TRUE],
                  sigma_d = u_d[, "mad", drop = TRUE],
                  mu_j = u_j[, "median", drop = TRUE],
                  sigma_j = u_j[, "mad", drop = TRUE])
  
  # Calculate model with Stan
  fit = stan(file = here("02_scripts", "04_counterfactual",
                         "03_counterfactual_model_general.stan"), 
             data = data_list, 
             iter = iter,
             chains = chains,
             seed = seed,
             verbose = TRUE)
  
  # Extract model parameter/weight distribution
  params = extract(fit)
  
  output = list(fit = fit,
                params = params)
  
  return(output)
}

### Functions for Weight Plotting
#### Function for correct naming of protected variables
param_names = function(col){
  
  pattern = c("_prot.1", 
              "_prot.2",
              "_prot.3",
              "_prot.4",
              "_prot.5",
              "_prot.6")
  
  replacement = c("_Race: Caucasian",
                  "_Race: African-American",
                  "_Race: Hispanic",
                  "_Race: Other",
                  "_Sex: Male",
                  "_Sex: Female")
  
  for (i in 1:length(pattern)){
    col = str_replace(col, pattern[i], replacement[i])
  }
  
  return(col)
}

#### Function for plot of compas weights (except age variable)
plot_compas_weights = function(data, plot_text){
  
  colors_discrete = rep(wes_palette("Darjeeling1"), 3)
  
  plot = data %>% 
    as.data.frame() %>% 
    select(starts_with("compas"),
           -c("compas_age")) %>% 
    rename_with(.cols = everything(), 
                .fn = ~ param_names(.x)) %>%
    rename_with(.cols = everything(), 
                .fn = ~ str_replace(.x, "compas_", "")) %>% 
    rename_with(.cols = everything(), 
                .fn = ~ str_replace(.x, "bias", "Bias")) %>% 
    rename_with(.cols = everything(), 
                .fn = ~ str_replace(.x, "ud", "Ud")) %>%   
    pivot_longer(cols = everything()) %>%
    ggplot(aes(x = name, y = value, fill = name)) +
    ggdist::stat_halfeye(
      adjust = 0.75,
      justification = -.05,
      .width = 0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = .12,
      outlier.colour = NA,
      alpha = 0.5
    ) +
    geom_abline(intercept = 0, slope = 0) +
    scale_fill_manual(values = colors_discrete) +
    guides(fill = "none") +
    labs(title = "Posterior-Verteilung der Modell-Parameter",
         subtitle = paste0("Nur Gewichte für COMPAS-Gleichung (ohne Alter) für ", plot_text),
         y = "Parameter-Wert",
         x = "Parameter für COMPAS-Prognose") +
    coord_flip() 
  
  return(plot)
}

#### Function for plot of compas weights (only age variable)
plot_compas_weights_age = function(data, plot_text, ymin = -0.5, ymax = 0){
  
  colors_discrete = rep(wes_palette("Darjeeling1"), 3)
  
  plot = data %>% 
    as.data.frame() %>% 
    select(c("compas_age")) %>% 
    rename_with(.cols = everything(), 
                .fn = ~ str_replace(.x, "compas_age", "Age")) %>% 
    rename_with(.cols = everything(), 
                .fn = ~ param_names(.x)) %>%
    rename_with(.cols = everything(), 
                .fn = ~ str_replace(.x, "compas_", "")) %>% 
    rename_with(.cols = everything(), 
                .fn = ~ str_replace(.x, "bias", "Bias")) %>% 
    rename_with(.cols = everything(), 
                .fn = ~ str_replace(.x, "ud", "Ud")) %>%   
    pivot_longer(cols = everything()) %>%
    ggplot(aes(x = name, y = value, fill = name)) +
    ggdist::stat_halfeye(
      adjust = 0.75,
      justification = -.05,
      .width = 0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = .12,
      outlier.colour = NA,
      alpha = 0.5
    ) +
    geom_abline(intercept = 0, slope = 0) +
    scale_fill_manual(values = colors_discrete) +
    guides(fill = "none") +
    ylim(ymin, ymax) +
    labs(title = "Posterior-Verteilung der Modell-Parameter",
         subtitle = paste0("Nur Gewichte für COMPAS-Gleichung (nur Alter) für ", plot_text),
         y = "Parameter-Wert",
         x = "Parameter für COMPAS-Prognose") +
    coord_flip() 
  
  return(plot)
}

#### Function for table of the median and mad of all weights
weight_table = function(data){
  data %>% 
    as.data.frame() %>% 
    select(-starts_with("u_d"), 
           -starts_with("u_j"), 
           -lp__) %>% 
    rename_with(.cols = everything(), 
                .fn = ~ param_names(.x)) %>% 
    summarise(across(everything(),
                     list(median = median,
                          mad = mad))) %>% 
    round(digits = 2) %>% 
    pivot_longer(cols = everything(),
                 names_to = c("variable",
                              "metric"),
                 names_pattern = "(.*)_([a-z]*)$") %>% 
    pivot_wider(names_from = "variable") %>% 
    column_to_rownames(var = "metric")
}

### Functions for Prediction & Plotting 
### (adapted code from https://medium.com/\@alex.pavlakis/making-predictions-from-stan-models-in-r-3e349dfac1ed

#### Function for prediction of logit values of given data frame
prediction_counter_logit = function(data, params_df, fun, seed){
  
  set.seed((seed))
  
  # Only U_d, U_j and compas weights are relevant
  params_df_new = params_df %>% 
    select(starts_with(c("u_", "compas_")))
  
  # Sample a parameter set for each observation
  params = data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(params_df_new)))
  colnames(params) = colnames(params_df_new)
  
  for (i in 1:ncol(params)) {
    params[,i] = sample(params_df_new[,i], size = nrow(params), replace = TRUE)
  }
  
  # Define certain parameter and variable data frames
  param_prot = params %>% 
    select(starts_with("compas_prot"))
  
  data_prot = data %>% 
    select(starts_with(c("Race_", "Sex_"))) %>% 
    as.matrix()
  
  data_ud = params %>% 
    select(starts_with("u_d")) %>% 
    as.matrix() %>% 
    t()
  
  # Create empty matrix for results
  inv_logit_raw = matrix(NA,
                         nrow = nrow(data),
                         ncol = ncol(data_ud)) 
  
  # Loop over all possible U_d values (4000 per observation)
  for (i in 1:ncol(data_ud)) {
    
    # Calculate linear combination for COMPAS equation for all observations
    lin_comb = params$compas_bias + 
      params$compas_age * data$Age + 
      rowSums(param_prot * data_prot) +
      as.vector(params$compas_ud * data_ud[,i])
    
    # Calculate logit for all observations
    inv_logit_raw[,i] = 1/(1 + exp(-lin_comb))
  }
  
  # Summarise logits over all possible U_d values with function "fun" (e.g. median)
  inv_logit = apply(inv_logit_raw, 1, fun)
  
  return(inv_logit)
}

#### Function for prediction dataframes and prediction plots (incl. weights) 
prediction_counter = function(params_df, 
                              data, 
                              threshold = 0.5, 
                              seed = 3928567, 
                              ymin_age_weight = -0.5,
                              ymax_age_weight = 0,
                              color_nr,
                              plot_text,
                              plot_text_bin,
                              plot_org_logit = NULL,
                              plot_org_bin = NULL){
  
  # Calculate predicted values for logits and COMPAS score
  counter_pred_logit = prediction_counter_logit(params_df = params_df %>% as.data.frame(),
                                                data = data,
                                                fun = "median",
                                                seed = seed)
  
  counter_pred = ifelse(counter_pred_logit >= threshold, 1, 0)
  
  # Calculate accuracy of prediction
  accuracy = mean((counter_pred == data$COMPAS))
  
  # Plots for weights
  plot_compas_weights = plot_compas_weights(data = params_df, 
                                            plot_text = plot_text)
  
  plot_compas_weights_age = plot_compas_weights_age(data = params_df, 
                                                    plot_text = plot_text, 
                                                    ymin = ymin_age_weight, 
                                                    ymax = ymax_age_weight)
  weight_table = weight_table(data = params_df)
  
  
  ## Plots for predictions
  colors = wes_palette("Darjeeling1")
  
  # Plot for predicted logits
  plot_counter_logit = counter_pred_logit %>% 
    as_tibble() %>% 
    rename(logit = value) %>% 
    ggplot(aes()) +
    geom_histogram(aes(x = logit, y = after_stat(density)), 
                   bins = 30, 
                   fill = colors[color_nr],
                   alpha = 0.3) +
    geom_density(aes(x = logit),
                 linewidth = 1,
                 color = colors[color_nr]) +
    labs(title = "Geschätzte Wahrscheinlichkeiten für einen niedrigen COMPAS-Score",
         subtitle = paste0("Histogramm und Dichtefunktion für ", plot_text),
         x = "Geschätzte Wahrscheinlichkeit für einen niedrigen COMPAS-Score",
         y = "Dichte") 
  
  # Plot for predicted COMPAS score
  plot_counter_bin = counter_pred %>% 
    as_tibble() %>% 
    group_by(value) %>% 
    summarise(n = n()) %>% 
    mutate(value = factor(value, levels = c(1,0), labels = c("Niedrig", "Hoch"))) %>% 
    ggplot(aes(x = value, y = n)) +
    geom_col(fill = colors[color_nr]) +
    labs(title = "Geschätzter Wert für den COMPAS-Score",
         subtitle = paste0("Absolute Häufigkeit der geschätzten Werte für ", plot_text),
         x = "COMPAS-Score",
         y = "Geschätzte Häufigkeit") 
  
  ## Plots for comparison of model predictions and base model predictions
  plot_comp_logit = NULL
  plot_comp_bin = NULL
  
  # Plot for comparison of predicted logits in model and base model
  if(!is.null(plot_org_logit)){
    plot_comp_logit = plot_org_logit +
      geom_histogram(data = counter_pred_logit %>% as_tibble() %>% rename(logit = value) ,
                     aes(x = logit, y = after_stat(density)), 
                     bins = 30, 
                     fill = colors[color_nr],
                     alpha = 0.3) +
      geom_density(data = counter_pred_logit %>% as_tibble() %>% rename(logit = value), 
                   aes(x = logit),
                   linewidth = 1,
                   color = colors[color_nr]) +
      labs(title = "Vergleich geschätzter Wahrscheinlichkeiten für einen niedrigen COMPAS-Score",
           subtitle = paste0("Histogramm und Dichtefunktion für originale & kontrafaktische Daten (", plot_text,")"),
           x = "Geschätzte Wahrscheinlichkeit für einen niedrigen COMPAS-Score",
           y = "Dichte") 
  }
  
  # Plot for comparison of predicted COMPAS scores in model and base model
  if(!is.null(plot_org_bin)){
    
    counter_pred_plot = counter_pred %>% 
      as_tibble() %>% 
      group_by(value) %>% 
      summarise(n = n()) %>% 
      mutate(Data = factor(plot_text_bin)) %>% 
      mutate(value = factor(value, levels = c(1,0), labels = c("Niedrig", "Hoch")))
    
    if(ncol(plot_org_bin$data) == 2){
      plot_comp_bin = plot_org_bin$data %>% 
        mutate(Data = factor("Original")) %>% 
        rbind(counter_pred_plot) %>% 
        group_by(Data) %>% 
        ggplot(aes(x = value, y = n, fill = Data)) +
        geom_col(position = position_dodge(0.82), width = 0.8) +
        labs(title = "Vergleich der geschätzten Werte für den COMPAS-Score",
             subtitle = "Häufigkeit der geschätzten Score-Werte für originale & kontrafaktische Daten",
             x = "COMPAS-Score",
             y = "Geschätzte Häufigkeit",
             fill = "Daten") +
        scale_fill_manual(values = c(colors[2], colors[color_nr])) 
    } else {
      plot_comp_bin = plot_org_bin$data %>% 
        rbind(counter_pred_plot) %>% 
        group_by(Data) %>% 
        ggplot(aes(x = value, y = n, fill = Data)) +
        geom_col(position = position_dodge(0.82), width = 0.8) +
        labs(title = "Vergleich der geschätzten Werte für den COMPAS-Score",
             subtitle = "Häufigkeit der geschätzten Score-Werte für originale & kontrafaktische Daten",
             x = "COMPAS-Score",
             y = "Geschätzte Häufigkeit",
             fill = "Daten") +
        scale_fill_manual(values = c(as.vector(unique(ggplot_build(plot_org_bin)$data[[1]]["fill"])[["fill"]]),
                                     colors[color_nr]))
    }
  }
  
  # Output list
  output = list(counter_pred_logit = counter_pred_logit,
                counter_pred = counter_pred,
                accuracy = accuracy,
                plot_compas_weights = plot_compas_weights,
                plot_compas_weights_age = plot_compas_weights_age,
                weight_table = weight_table,
                plot_counter_logit = plot_counter_logit,
                plot_counter_bin = plot_counter_bin,
                plot_comp_logit = plot_comp_logit,
                plot_comp_bin = plot_comp_bin
  )
  
  return(output)
}

### Function for approximate counterfactual fairness
approx_counterfairness = function(base_outcome = pred_base$counter_pred,
                                  counter_outcome,
                                  epsilon){
  
  diff = abs(base_outcome-counter_outcome)
  
  p = diff < epsilon
  
  delta = 1 - sum(p)/length(p)
  
  return(delta)
}

## Functions for Post-Processing
### Functions for outcome calculation
normal_function = function(logit){
  output = ifelse(logit >= 0.5, 1, 0)
  return(output)
}

rejected_function = function(data, privileged){
  
  d = data %>% 
    select(all_of(privileged))
  
  output = ifelse(d == 1, 0, 1)
  
  return(output)
}

### Delta calculation given a fixed theta and epsilon
delta_calc = function(theta,
                      epsilon = 0.1,
                      data_base = data_base,
                      logit_base,
                      data_counter, 
                      logit_counter, 
                      prot, 
                      priv){
  
  upper = 0.5 + theta
  lower = 0.5 - theta
  
  d_base = data_base %>% 
    mutate(logit = logit_base) %>% 
    select(starts_with(prot),
           logit)
  
  d_base_outcome = d_base %>% 
    mutate(outcome = ifelse((logit < lower | logit > upper),
                            normal_function(logit),
                            rejected_function(data = ., privileged = priv))
    ) 
  
  deltas = c()
  
  for (i in 1:length(data_counter)) {
    d_counter = data_counter[[i]] %>% 
      mutate(logit = logit_counter[[i]]) %>% 
      select(starts_with(prot),
             logit)
    
    d_counter_outcome = d_counter %>%
      mutate(outcome = ifelse((logit < lower | logit > upper),
                      normal_function(logit),
                      rejected_function(data = ., privileged = priv))
             )
    
    # d_counter_outcome = d_counter %>%
    #   mutate(outcome = normal_function(logit))
    
    deltas[i] = approx_counterfairness(base_outcome = d_base_outcome$outcome,
                                       counter_outcome = d_counter_outcome$outcome,
                                       epsilon = epsilon)
  }
  
  delta = max(deltas)
  
  return(delta)
}

### Theta calculation given a fixed range of possible thetas
theta_calc = function(theta_lower = 0,
                      theta_upper = 0.05,
                      step = 0.01,
                      epsilon = 0.1,
                      data_base = data_base,
                      logit_base,
                      data_counter, 
                      logit_counter, 
                      prot, 
                      priv){
  
  grid = seq(theta_lower, theta_upper, step)
  
  df = data.frame(theta = grid,
                  delta = NA)
  
  for(i in 1:length(grid)){
    df$delta[i] = delta_calc(theta = df$theta[i],
                             epsilon = epsilon,
                             data_base = data_base,
                             logit_base = logit_base,
                             data_counter = data_counter, 
                             logit_counter = logit_counter, 
                             prot = prot, 
                             priv = priv)
  }
  
  output = df %>% 
    filter(delta == min(delta)) %>% 
    as.list()
  
  return(output)
}
