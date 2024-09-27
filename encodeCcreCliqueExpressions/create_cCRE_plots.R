library(ggplot2)
library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(reshape2)
library(ggstatsplot)

datapchic <- read.csv("results/stats-encodeCCRE-tisPCHiC.csv", header=T, dec=".")
datahic <- read.csv("results/stats-encodeCCRE-3DIV.csv", header=T, dec=".")

datapchic1 <- filter(datapchic,mean_FANTOM5_for_all_nodes>0,cardinality_C5segs_FANTOM5_for_all_segments_that_are_in_any_clique5>50)
datapchic2 <- filter(datapchic,mean_GTEx_for_all_nodes>0,cardinality_C5segs_GTEx_for_all_segments_that_are_in_any_clique5>50)
datahic1 <- filter(datahic,mean_FANTOM5_for_all_nodes>0,cardinality_C5segs_FANTOM5_for_all_segments_that_are_in_any_clique5>50)
datahic2 <- filter(datahic,mean_GTEx_for_all_nodes>0,cardinality_C5segs_GTEx_for_all_segments_that_are_in_any_clique5>50)

#Figure 1

datapchic_medgrank <- melt(datapchic2,id.vars=c('chr','tissue'),
                           measure.vars=c('med_GTEx_for_all_nodes',
                                          'med_C3segs_GTEx_for_all_segments_that_are_in_any_clique',
                                          'med_C4segs_GTEx_for_all_segments_that_are_in_any_clique4',
                                          'med_C5segs_GTEx_for_all_segments_that_are_in_any_clique5'))

ggbetweenstats(
  data = datapchic_medgrank,
  x = variable,
  y = value,
  type = "nonparametric",
  plot.type = "box",
  xlab = "Clique size",
  ylab = "Median GTEx expression",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  title = "Clique size effect on GTEx expression in pcHiC data",
  ggplot.component = scale_x_discrete(labels=c('Single node','C3 clique','C4 clique','C5 clique'))
)

datapchic_medfrank <- melt(datapchic1,id.vars=c('chr','tissue'),
                           measure.vars=c('med_FANTOM5_for_all_nodes',
                                          'med_C3segs_FANTOM5_for_all_segments_that_are_in_any_clique',
                                          'med_C4segs_FANTOM5_for_all_segments_that_are_in_any_clique4',
                                          'med_C5segs_FANTOM5_for_all_segments_that_are_in_any_clique5'))

ggbetweenstats(
  data = datapchic_medfrank,
  x = variable,
  y = value,
  type = "nonparametric",
  plot.type = "box",
  xlab = "Clique size",
  ylab = "Median FANTOM5 expression",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  title = "Clique size effect on FANTOM5 expression in pcHiC data",
  ggplot.component = scale_x_discrete(labels=c('Single node','C3 clique','C4 clique','C5 clique'))
)

#Figure 2

datapchic_medgrankE <- melt(datapchic2,id.vars=c('chr','tissue'),
                            measure.vars=c('med_E_GTEx_for_all_segments_with_distal_enhancers',
                                           'med_EE_links_GTEx_for_all_segments_that_are_in_any_link_where_both_nodes_have_these_encode_props',
                                           'EEE_med_GTEx',
                                           'med_EEEE_GTEx_for_all_segments_that_are_in_a_clique4_where_all_4_nodes_have_these_encode_props',
                                           'med_EEEEE_GTEx_for_all_segments_that_are_in_a_clique5_where_all_5_nodes_have_these_encode_props'))

datapchic_medgrankE1 <- filter(datapchic_medgrankE,value>-1)

ggbetweenstats(
  data = datapchic_medgrankE1,
  x = variable,
  y = value,
  type = "nonparametric",
  plot.type = "box",
  xlab = "Enhancer clique size",
  ylab = "Median GTEx expression",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  title = "Clique size effect on expression in enhancer nodes in pcHi-C",
  ggplot.component = scale_x_discrete(labels=c('E','EE','EEE','EEEE','EEEEE'))
)

datapchic_medgC3E <- melt(datapchic2,id.vars=c('chr','tissue'),
                          measure.vars=c('med_C3segs_GTEx_for_all_segments_that_are_in_any_clique',
                                         'med_C3EncS_GTEx_for_all_segments_that_are_in_a_clique_where_all_3_nodes_have_some_encode_prop',
                                         'med_C3_E_GTEx_clique_nodes_with_a_distal_enhancer_in_a_set_of_cliques_where_all_3_nodes_have_some_encode_prop',
                                         'EE._med_GTEx',
                                         'EEE_med_GTEx'))

datapchic_medgC3E1 <- filter(datapchic_medgC3E,value>-1)

ggbetweenstats(
  data = datapchic_medgC3E1,
  x = variable,
  y = value,
  type = "nonparametric",
  plot.type = "box",
  xlab = "Number of distal enhancer nodes in clique",
  ylab = "Median GTEx expression",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  title = "Distal enhancer concentration in pcHi-C cliques",
  ggplot.component = scale_x_discrete(labels=c('All C3','All cCRE','E**','EE*','EEE')))

#Figure 3

datapchic_medgC3PE <- melt(datapchic2,id.vars=c('chr','tissue'),
                         measure.vars=c('med_C3_P_GTEx_clique_nodes_with_a_promoter_in_a_set_of_cliques_where_all_3_nodes_have_some_encode_prop',
                                        'PPP_med_GTEx',
                                        'PPE_med_GTEx',
                                        'PEE_med_GTEx',
                                        'EEE_med_GTEx',
                                        'med_C3_E_GTEx_clique_nodes_with_a_distal_enhancer_in_a_set_of_cliques_where_all_3_nodes_have_some_encode_prop'))

datapchic_medgC3PE1 <- filter(datapchic_medgC3PE,value>-1)

ggbetweenstats(
  data = datapchic_medgC3PE1,
  x = variable,
  y = value,
  type = "nonparametric",
  plot.type = "box",
  xlab = "Clique cCRE composition",
  ylab = "Median GTEx expression",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  title = "Clique annotation combination effect on GTEx expression in pcHi-C",
  ggplot.component = scale_x_discrete(labels=c('P**','PPP','PPE','PEE','EEE','E**'))
)
