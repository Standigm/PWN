# Copyright 2022 Standigm Inc.. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

library(tidyverse)


#' Load files from various PPI.
#'
#' @param filename Filename without extension.
#' @param sources Sources of PPI.
#' @return Tibble
read_files <- function(filename, sources = c("STRING", "BioGRID")) {
  sources %>%
    map_dfr(~ read_tsv(str_c(filename, ".", str_to_lower(.x), ".tsv"), na = character()) %>% mutate(ppi = .x)) %>%
    mutate(ppi = factor(ppi, levels = sources))
}

#' Save the last figure to given location.
#'
#' @param filename Filename without extension.
#' @param height Relative height scale factors.
#' @param ... Other arguments passed on to \code{theme}.
#' @return None
plot_save <- function(filename, height = 1, singlewidth = FALSE, ...) {
  width <- ifelse(singlewidth, 86, 178)
  for (ext in c("pdf", "png")) {
    ggsave(str_c("figure/", filename, ".", ext),
      last_plot() + theme_bw() + theme(legend.position = "bottom", ...),
      width = width, height = height * width, unit = "mm", dpi = 350
    )
  }
}

#' Prettify p-value.
#'
#' @param pvalue p-value
#' @return Prettified p-value
format_pvalue <- function(pvalue) {
  pvalue %>%
    is.na() %>%
    if_else("---", formatC(pvalue, format = "e", digit = 2))
}

pvalue_str <- function(test_result) {
  pvalue <- test_result$p.value
  ifelse(pvalue > 0.0,
    str_c("= ", formatC(pvalue, format = "e", digits = 2)),
    str_c("< ", formatC(.Machine$double.eps, format = "e", digits = 2))
  )
}

set.seed(42)

rho0 <- 0.2
beta0 <- 0.5
gamma0 <- 0.5
method_keys <- c(
  "Naive",
  "RWR",
  "RWR w/ GDC",
  "mND",
  "RWR + Curv",
  "RWR + Prior(uKIN)",
  "RWR + Prior(uKIN, Curv)",
  "RWR + Prior(uKIN) + Curv",
  "RWR + Prior(uKIN, Curv) + Curv",
  "RWR + Prior(RWR)",
  "RWR + Prior(RWR, Curv)",
  "RWR + Prior(RWR) + Curv",
  "RWR + Prior(RWR, Curv) + Curv"
)
method_labels <- c(
  "Naive",
  "RWR",
  "RWR w/ GDC",
  "mND",
  "RWR w/ curvature",
  "uKIN",
  "RWR + Prior(uKIN, curvature)",
  "RWR + Prior(uKIN) + curvature",
  "uKIN w/ curvature",
  "PWN w/o curvature",
  "RWR + Prior(RWR, curvature)",
  "RWR + Prior(RWR) + curvature",
  "PWN"
)


weights <- read_files("edge") %>%
  mutate(beta = (round(100 * na_if(beta, -Inf)) / 100) %>% as_factor())
weights %>%
  ggplot() +
  geom_histogram(aes(x = weight, y = ..density.., fill = beta), bins = 50) +
  facet_grid(ppi ~ beta, scale = "free_y", labeller = labeller(.cols = NULL)) +
  scale_fill_brewer(palette = "RdYlBu") +
  labs(x = expression("Edge weight of" ~ K), fill = "beta")
plot_save("weights", height = 0.67, strip.text.x = element_blank())

topo <- read_files("topo") %>% mutate(degree = as.integer(degree))
topo_kde <- topo %>%
  group_by(ppi) %>%
  group_modify(function(df, key) {
    h <- 50
    n <- 100
    eps <- 1e-6
    x <- df %>% pull(degree)
    y <- df %>% pull(curv_mean)
    area <- c(min(x) - h / 2, max(x) + h / 2, min(y) - h / 2, max(y) + h / 2)
    density <- MASS::kde2d(x, y, h = h, n = n, lims = area)
    x_prior <- df %>%
      filter(prior) %>%
      pull(degree)
    y_prior <- df %>%
      filter(prior) %>%
      pull(curv_mean)
    density_prior <- MASS::kde2d(x_prior, y_prior, h = h, n = n, lims = area)
    tibble(
      size = max(area[2] - area[1], area[4] - area[3]) / n,
      degree = density_prior$x %>% rep(times = length(density_prior$y)),
      curv_mean = density_prior$y %>% rep(each = length(density_prior$x)),
      overall = c(density$z),
      prior = c(density_prior$z),
      relative = c((eps + density_prior$z) / (eps + density$z))
    ) %>%
      filter(overall > eps * 0.01)
  }) %>%
  gather(overall, prior, relative, key = "type", value = "value") %>%
  mutate(type = factor(
    type, c("overall", "prior", "relative"),
    c(
      "Density of all nodes", "Density of prior nodes",
      "Relative log density of prior"
    )
  )) %>%
  ungroup()
topo_kde %>%
  filter(type == "Density of prior nodes") %>%
  ggplot(aes(x = degree, y = curv_mean)) +
  geom_tile(aes(fill = value, height = size, width = size)) +
  scale_fill_distiller(
    palette = "RdYlBu", direction = -1, trans = "sqrt",
    breaks = c(1e-6, 2e-5, 5e-5)
  ) +
  guides(fill = guide_colourbar(barwidth = 12)) +
  coord_fixed(ratio = 1) +
  facet_wrap(~ppi) +
  labs(
    x = "Node degree", y = "Average curvature", fill = "Relative density of priors"
  )
plot_save("eda", height = 0.67)

var <- read_files("var")
var %>%
  group_by(ppi) %>%
  group_map(~ broom::tidy(lm(uKIN ~ RWR - 1, data = .x))) %>%
  print()
var %>%
  ggplot(aes(x = RWR, y = uKIN)) +
  stat_bin_2d(aes(fill = stat(density)), bins = 50) +
  geom_hline(aes(yintercept = 0),
    color = RColorBrewer::brewer.pal(5, "RdGy")[4], size = 0.25
  ) +
  geom_vline(aes(xintercept = 0),
    color = RColorBrewer::brewer.pal(5, "RdGy")[4], size = 0.25
  ) +
  geom_abline(aes(intercept = 0, slope = 1),
    color = RColorBrewer::brewer.pal(5, "RdGy")[5], size = 0.25, linetype = "42"
  ) +
  facet_wrap(~ppi) +
  scale_fill_distiller(trans = "log10", palette = "RdYlBu") +
  coord_fixed(ratio = 1.0) +
  labs(
    x = "Standard deviation of PWN", y = "Standard deviation of uKIN",
    color = "", fill = "Density"
  ) +
  guides(x = guide_axis(angle = 45), y = guide_axis(angle = 45), fill = guide_colourbar(barwidth = 10))
plot_save("exp5", height = 0.75)


rawdf <- read_files("perf") %>%
  mutate(
    method = factor(method, method_keys, labels = method_labels),
    beta = round(100 * na_if(beta, -Inf)) / 100,
    gamma = round(100 * na_if(gamma, -Inf)) / 100,
    rho = round(100 * na_if(rho, -Inf)) / 100
  )


df1 <- rawdf %>%
  filter(
    method %in% c(
      "RWR",
      "RWR w/ GDC",
      "mND",
      "uKIN",
      "RWR w/ curvature",
      "uKIN w/ curvature",
      "PWN w/o curvature",
      "PWN"
    ),
    rho == rho0,
    (is.na(beta) | beta == beta0),
    (is.na(gamma) | gamma == gamma0)
  ) %>%
  filter(method != "Naive") %>%
  select(ppi, idx, method, metric, perf)

df1 %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = perf, fill = method), outlier.size = 0.5) +
  facet_wrap(~ ppi + metric, scale = "free_y") +
  scale_fill_brewer(palette = "RdYlBu", na.translate = FALSE) +
  labs(x = "", y = "Performance", fill = "Method") +
  guides(fill = guide_legend(nrow = 4))
plot_save("exp1", axis.text.x = element_blank())

baseline <- df1 %>%
  filter(method == "RWR") %>%
  transmute(ppi, idx, metric, baseline = perf)
df1 %>%
  inner_join(baseline) %>%
  group_by(ppi, method, metric) %>%
  summarise(
    mean = mean(perf), sd = sd(perf),
    pvalue = t.test(perf - baseline,
      alternative = "greater"
    )$p.value
  ) %>%
  ungroup() %>%
  transmute(ppi, method, metric,
    result = sprintf("%f (%s)", mean, format_pvalue(pvalue))
  ) %>%
  spread(metric, result) %>%
  write_csv(str_c("table/exp1.csv"))


df2 <- rawdf %>%
  filter(is.na(rho)) %>%
  filter(method == "RWR w/ curvature")

df2 %>%
  ggplot() +
  geom_hline(aes(yintercept = perf, color = "RWR"),
    data = df2 %>% filter(beta == 0.0)
  ) +
  geom_line(aes(x = beta, y = perf, color = "RWR w/ curvature")) +
  scale_color_manual(
    values = c(
      "RWR" = RColorBrewer::brewer.pal(4, "RdYlBu")[1],
      "RWR w/ curvature" = RColorBrewer::brewer.pal(4, "RdYlBu")[4]
    ),
    breaks = c("RWR", "RWR w/ curvature")
  ) +
  facet_wrap(~ ppi + metric) +
  labs(x = "Beta", y = "Performance", color = "Method")
plot_save("exp2")

df3 <- rawdf %>%
  filter(rho == rho0, method %in% c(
    "PWN",
    "uKIN"
  )) %>%
  mutate(
    beta = beta %>%
      as_factor() %>%
      fct_relabel(~ if_else(.x == 0, "No curv", sprintf("beta = %s", .x))),
    gamma = gamma %>% as_factor() %>% fct_explicit_na("uKIN")
  )

df3 %>%
  ggplot() +
  geom_boxplot(aes(x = gamma, y = perf, fill = beta), outlier.size = 0.5) +
  facet_wrap(~ ppi + metric, scale = "free_y", nrow = 3, dir = "v") +
  scale_fill_brewer(palette = "RdYlBu", na.translate = FALSE) +
  labs(x = "Restart probabilty", y = "Performance", fill = "PWN")
plot_save("exp3", height = 1.25)


df4 <- rawdf %>%
  filter(!is.na(rho)) %>%
  filter(method == "RWR" |
    method == "uKIN" |
    (method == "RWR w/ curvature" & beta == beta0) |
    (method == "PWN" &
      beta == beta0 & gamma == gamma0)) %>%
  mutate(
    rho = if_else(method == "RWR" | method == "RWR w/ curvature",
      0, 100 * rho
    ) %>% as_factor(),
    method = method %>%
      fct_recode(
        `uKIN` = "RWR",
        `PWN` = "RWR w/ curvature"
      )
  )

df4 %>%
  ggplot() +
  geom_boxplot(aes(x = rho, y = perf, fill = method), outlier.size = 0.5) +
  facet_wrap(~ ppi + metric, scale = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "% of train set", y = "Performance", fill = "Method")
plot_save("exp4")
