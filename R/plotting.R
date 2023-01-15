### ZICRM-ASCA+ plotting
### Author: Auke Haver
### BDA GROUP SILS Amsterdam
### Project: Zero-inflated GLMM ASCA

plot_simulated_time_scores <- function(design_matrix, score_matrix, expl_var, component, y_label = ""){
  ggplot() +
    geom_point(
      data = design_matrix %>%
        mutate(score = as.matrix(design_matrix[,5:6]) %*% Ta[,component]) %>%
        distinct(timepoint,score),
      mapping = aes(x = timepoint,y = score),
      shape = 21,
      color = "black",
      fill = "white"
    )+
    # Theme
    theme_bw()+
    #theme_classic(base_family="Times")+
    # Labels
    labs(
      x = "Timepoint",
      y=substitute(paste(m, " PC", c, " (",e,"%)", sep = ""), list(m = y_label, c = component, e = expl_var[component])),
    )
}
plot_simulated_interaction_scores <- function(design_matrix, score_matrix, expl_var, component, y_label = ""){
  ggplot(
    data = design_matrix %>%
      mutate(score = as.matrix(design_matrix[,7:14]) %*% Tab[,component]) %>%
      distinct(timepoint, treatment, score)
    ) +
    geom_point(
      mapping = aes(x = timepoint, y = score, fill = factor(treatment)),
      shape = 21
    )+
    geom_line(
      # data = design_matrix %>%
      #   mutate(score = score_matrix[,component]) %>%
      #   distinct(timepoint, treatment, score),
      mapping = aes(
        x = timepoint,
        y = score,
        group = factor(treatment),
        color = factor(treatment)
      ),
      show.legend = FALSE
    ) +
    # Theme and labels
    theme_bw()+
    #theme(axis.text.x=element_text(size=10))+
    scale_color_viridis_d()+
    scale_fill_viridis_d()+
    labs(
      y=substitute(paste(m, " PC", c, " (",e,"%)", sep = ""), list(m = y_label, c = component, e = expl_var[component])),
      x="Timepoint",
      fill="Treatment Group"
    )
}
plot_simulated_loadings <- function(loading_matrix, expl_var, component, y_label = ""){
  ggplot() +
    geom_col(
      data = tibble(
        var = factor(seq(1:nrow(loading_matrix))),
        loading = loading_matrix[,component]
      ),
      mapping = aes(
        x = LETTERS[var],
        y = loading
      ),
      color = "black",
      fill = "white"
    )+
    theme_bw()+
    # Theme and labels
    #theme_classic(base_family = "Times")+
    #theme(axis.text.x=element_text(size=6))+
    labs(
      y=substitute(paste(m, " PC", c, " (",e,"%)", sep = ""), list(m = y_label, c = component, e = expl_var[component])),
      x="Variable"
    )
}
plot_estimate_time_scores <- function(score_dataframe, perc_expl_dataframe, component, y_label = ""){
  expl_var <- perc_expl_dataframe %>%
    filter(
      fold != "ref",
      PC == component
    ) %>%
    pull("perc_expl")
  ggplot() +
    geom_errorbar(
      data = score_dataframe %>%
        filter(
          fold != "ref",
          PC == component
        ) %>%
        group_by(timepoint) %>%
        summarize(
          min_score = min(score),
          max_score = max(score),
          .groups = "drop"
        ),
      mapping = aes(
        x = timepoint,
        ymin = min_score,
        ymax = max_score
      ),
      show.legend = FALSE
    )+
    geom_point(
      data = score_dataframe %>%
        filter(
          fold == "ref",
          PC == component
        ),
      mapping = aes(
        x = timepoint,
        y = score
      )
    ) +
    theme_bw()+
    # Theme and labels
    #theme_classic(base_family = "Times")+
    #theme(axis.text.x=element_text(size=10))+
    scale_color_viridis_d()+
    scale_fill_viridis_d()+
    labs(
      y = substitute(paste(m, " PC", c, " (", min, "-", max, "%)", sep = ""), list(m = y_label, min = format(min(expl_var),nsmall=1), max = format(max(expl_var),nsmall=1), c = component)),
      #y=            paste0("Scores PC", component, ),
      x="Timepoint"
    )
}
plot_estimate_interaction_scores <- function(score_dataframe, perc_expl_dataframe, component, y_label = ""){
  expl_var <- perc_expl_dataframe %>%
    filter(
      fold != "ref",
      PC == component
    ) %>%
    pull("perc_expl")
  ggplot() +
    geom_line(
      data = score_dataframe %>%
        filter(
          fold == "ref",
          PC == component
        ),
      mapping = aes(
        x = timepoint,
        y = score,
        color = factor(treatment),
        group = factor(treatment)
      ),
      show.legend = FALSE
    )+
    geom_errorbar(
      data = score_dataframe %>%
        filter(
          fold != "ref",
          PC == component
        ) %>%
        group_by(timepoint, treatment) %>%
        summarize(
          min_score = min(score),
          max_score = max(score),
          .groups = "drop"
        ),
      mapping = aes(
        x = timepoint,
        ymin = min_score,
        ymax = max_score,
        color = factor(treatment)
      ),
      show.legend = FALSE
    )+
    geom_point(
      data = score_dataframe %>%
        filter(
          fold == "ref",
          PC == component
        ),
      mapping = aes(
        x = timepoint,
        y = score,
        fill = factor(treatment),
        group = factor(treatment)
      ),
      
      shape = 21
    ) +
    theme_bw()+
    # Theme and labels
    #theme_classic(base_family = "Times")+
    #theme(axis.text.x=element_text(size=10))+
    scale_color_viridis_d()+
    scale_fill_viridis_d()+
    labs(
      y = substitute(paste(m, " PC", c, " (", min, "-", max, "%)", sep = ""), list(m = y_label, min = format(min(expl_var),nsmall=1), max = format(max(expl_var),nsmall=1), c = component)),
      #y=paste0("Scores PC", component, " (", min(perc_expl_list), "-", max(perc_expl_list), "%)"),
      x="Timepoint",
      fill="Treatment Group"
    )
}
plot_estimate_loadings <- function(loadings_dataframe, perc_expl_dataframe, component, y_label = ""){
  expl_var <- perc_expl_dataframe %>%
    filter(
      fold != "ref",
      PC == component
    ) %>%
    pull("perc_expl")
  ggplot() +
    geom_col(
      data = loadings_dataframe %>%
        filter(
          PC == component,
          fold == "ref"
        ),
      mapping = aes(
        x = LETTERS[as.numeric(var)],
        y = loading
      ),
      color = "black",
      fill = "white"
    ) +
    geom_errorbar(
      data = loadings_dataframe %>%
        filter(
          PC == component,
          fold != "ref"
        ) %>%
        group_by(var) %>%
        summarize(
          min_loading = min(loading),
          max_loading = max(loading),
          .groups = "drop"
        ),
      mapping = aes(
        x = LETTERS[var],
        ymin = min_loading,
        ymax = max_loading
      )
    )+
    geom_hline(
      mapping = aes(
        yintercept = 0
      ),
      linetype = "dashed"
    ) +
    theme_bw()+
    # Theme and labels
    #theme_classic(base_family = "Times")+
    #theme(axis.text.x=element_text(size=6))+
    labs(
      y = substitute(paste(m, " PC", c, " (", min, "-", max, "%)", sep = ""), list(m = y_label, min = format(min(expl_var),nsmall=1), max = format(max(expl_var),nsmall=1), c = component)),
      #y=paste0("Loadings PC", component, " (", min(perc_expl_list), "-", max(perc_expl_list), "%)"),
      x="Variable"
    )
}
