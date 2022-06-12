### ZICRM-ASCA+ plotting
### Author: Auke Haver
### BDA GROUP SILS Amsterdam
### Project: Zero-inflated GLMM ASCA

plot_simulated_time_scores <- function(design_matrix, score_matrix, expl_var, component){
  ggplot() +
    geom_point(
      data = design_matrix %>%
        mutate(score = score_matrix[,component]) %>%
        distinct(
          timepoint,
          score
        ),
      mapping = aes(
        x = timepoint,
        y = score
      ),
      shape = 21,
      color = "black",
      fill = "white"
    )+
    # Theme
    theme_classic(base_family="Times")+
    # Labels
    labs(
      x = "Timepoint",
      y = paste0(
        "Scores PC", 
        component, 
        " (",
        expl_var[component],
        "%)"
      )
    )
}
plot_simulated_interaction_scores <- function(design_matrix, score_matrix, expl_var, component){
  ggplot() +
    geom_point(
      data = design_matrix %>%
        mutate(score = score_matrix[,component]) %>%
        distinct(
          timepoint,
          treatment,
          score
        ),
      mapping = aes(
        x = timepoint,
        y = score,
        fill = treatment
      ),
      shape = 21
    )+
    geom_line(
      data = design_matrix %>%
        mutate(score = score_matrix[,component]) %>%
        distinct(
          timepoint,
          treatment,
          score
        ),
      mapping = aes(
        x = timepoint,
        y = score,
        group = treatment,
        color = treatment
      ),
      show.legend = FALSE
    ) +
    # Theme and labels
    theme_classic(base_family = "Times")+
    theme(axis.text.x=element_text(size=10))+
    scale_color_viridis_d()+
    scale_fill_viridis_d()+
    labs(
      y=paste0("Scores PC", component, " (", expl_var[component], "%)"),
      x="Timepoint",
      fill="Treatment Group"
    )
}
plot_simulated_loadings <- function(loading_matrix, expl_var, component){
  ggplot() +
    geom_col(
      data = tibble(
        var = factor(seq(1:nrow(loading_matrix))),
        loading = loading_matrix[,component]
      ),
      mapping = aes(
        x = reorder(var, as.numeric(var)),
        y = loading
      ),
      color = "black",
      fill = "white"
    )+
    # Theme and labels
    theme_classic(base_family = "Times")+
    theme(axis.text.x=element_text(size=6))+
    labs(
      y=paste0("Loadings PC", component, " (", expl_var[component], "%)"),
      x="Variable"
    )
}
plot_estimate_time_scores <- function(score_dataframe, perc_expl_dataframe, component){
  perc_expl_list <- perc_expl_dataframe %>%
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
    # Theme and labels
    theme_classic(base_family = "Times")+
    theme(axis.text.x=element_text(size=10))+
    scale_color_viridis_d()+
    scale_fill_viridis_d()+
    labs(
      y=paste0("Scores PC", component, " (", min(perc_expl_list), "-", max(perc_expl_list), "%)"),
      x="Timepoint"
    )
}
plot_estimate_interaction_scores <- function(score_dataframe, perc_expl_dataframe, component){
  perc_expl_list <- perc_expl_dataframe %>%
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
        color = treatment,
        group = treatment
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
        color = treatment
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
        fill = treatment,
        group = treatment
      ),
      
      shape = 21
    ) +
    # Theme and labels
    theme_classic(base_family = "Times")+
    theme(axis.text.x=element_text(size=10))+
    scale_color_viridis_d()+
    scale_fill_viridis_d()+
    labs(
      y=paste0("Scores PC", component, " (", min(perc_expl_list), "-", max(perc_expl_list), "%)"),
      x="Timepoint",
      fill="Treatment Group"
    )
}
plot_estimate_time_loadings <- function(loadings_dataframe, perc_expl_dataframe, component){
  perc_expl_list <- perc_expl_dataframe %>%
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
        x = reorder(var, as.numeric(var)),
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
        x = var,
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
    # Theme and labels
    theme_classic(base_family = "Times")+
    theme(axis.text.x=element_text(size=6))+
    labs(
      y=paste0("Loadings PC", component, " (", min(perc_expl_list), "-", max(perc_expl_list), "%)"),
      x="Variable"
    )
}