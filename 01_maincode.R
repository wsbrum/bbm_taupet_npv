library(tidyverse)
library(caret)
library(OptimalCutpoints)
library(glue)
library(aba)
library(ggdist)
library(ggh4x)


dementia_dx <-
  c(
    'AD',
    "Dementia_bvFTD",
    "Dementia_DLB",
    "Dementia_FTD_NOS",
    "Dementia_Not_determined",
    "Dementia_svPPA",
    "Dementia_VaD",
    "Dementia_Neurodegenerative_disorder_NOS",
    "Dementia_PD",
    "Dementia_PDD",
    "Dementia_nfvPPA",
    "Dementia_non_neurodegenerative",
    "Dementia_Parkinsonism_NOS",
    "Dementia_CBS"
  )

##################
# Basic handling #
##################

data <- readr::read_csv('data/df_processed.csv') %>%
  mutate(plasma_ab4240_araclon = -1*plasma_ab4240_araclon) %>% # For equivalent-direction interpretability
  select(
    sid,
    dx_baseline, age, sex_baseline,
    TAU_PET_STATUS,
    plasma_ab4240_araclon,
    plasma_ptau181_lilly,
    plasma_ptau217_lilly,
    plasma_ptau231_ugot,
    plasma_nfl_simoa,
    plasma_gfap_simoa
  ) %>%
  filter(
    !is.na(dx_baseline), !is.na(age),
    !is.na(sex_baseline),
    !is.na(TAU_PET_STATUS),
    dx_baseline %in% c('SCD','MCI', dementia_dx)
  )

data <- data %>%
  group_by(sid) %>%
  mutate(
    has_plasma = !all(
      is.na(plasma_ab4240_araclon) &
        is.na(plasma_ptau181_lilly) &
        is.na(plasma_ptau217_lilly) &
        is.na(plasma_ptau231_ugot) &
        is.na(plasma_nfl_simoa) &
        is.na(plasma_gfap_simoa)
    )
  ) %>%
  ungroup() %>%
  filter(has_plasma)

####################################################################################
# Converting concentrations to probabilities for better cross-assay interpretation #
####################################################################################

model <- data %>%
  aba_model() %>%
  set_groups(
    dx_baseline %in% c('SCD','MCI', dementia_dx),
    dx_baseline %in% c('SCD','MCI'),
    dx_baseline %in% dementia_dx,
    .labels = c('All participants', 'SCD-MCI', 'All-cause dementia')
  ) %>%
  set_outcomes(
    TAU_PET_STATUS
  ) %>%
  # set_covariates(
  #   age, gender_baseline_variable,
  #   .include_basic=F
  # ) %>%
  set_predictors(
    plasma_ab4240_araclon,
    plasma_ptau181_lilly,
    plasma_ptau217_lilly,
    plasma_ptau231_ugot,
    plasma_nfl_simoa,
    plasma_gfap_simoa,
    .labels = c(
      'Ab42/Ab40', 'pTau181', 'pTau217', 'pTau231', 'NfL', 'GFAP'
    )
  ) %>%
  set_stats(
    stat_glm(std.beta=T, complete.cases=F)
  ) %>%
  set_evals(
    eval_boot(ntrials=1000)
  ) %>%
  fit()


##########################
# Get demographics table #
##########################

model %>%
  aba_demographics(strata='dx_baseline') %>%
  aba_write('manuscript/Table1.csv')


################################
# Get individual probabilities #
################################

pred_df <- model %>% aba_predict(augment=T, merge=F)

sens_cut_fun <- function(sens_val, data) {

  sens_val <- sens_val
  cut_val <- quantile(filter(data, TAU_PET_STATUS ==1)$.fitted,
                      1 - sens_val,
                      type = 4)

  Predicted <- factor(
    as.integer(
      data$.fitted > cut_val
    ),
    levels=c(1,0),
    labels=c('Disease','Healthy')
  )
  Truth <- factor(as.integer(data$TAU_PET_STATUS),
                  levels=c(1,0),
                  labels=c('Disease','Healthy'))
  c <- table(Predicted, Truth)

  cc <- caret::confusionMatrix(c)
  ppv <- cc$byClass[['Pos Pred Value']]
  npv <- cc$byClass[['Neg Pred Value']]
  se <- cc$byClass[['Sensitivity']]
  sp <- cc$byClass[['Specificity']]
  prev <- cc$byClass[['Prevalence']]
  percent_under <- mean(Predicted == 'Healthy')

  data.frame(
    sens_val = sens_val,
    cut_val = cut_val,
    prev = prev,
    npv = npv,
    ppv = ppv,
    se = se,
    sp = sp,
    pu = percent_under
  )
}

# Set a range of sensitivity values
sens_vals <- c(0.5,0.6,0.7,0.8,0.9,0.95,0.975,0.99)

res_df <- pred_df %>%
  rowwise() %>%
  mutate(
    res = list(sens_vals %>% map_df(sens_cut_fun, data=.data$data))
  ) %>%
  unnest(res) %>%
  select(-data)
res_df0 <- res_df
res_df <- res_df0 %>%
  filter(sens_val < se)

# sds
res_sd <- res_df %>%
  filter(trial != 0) %>%
  group_by(group, outcome, stat, predictor, sens_val) %>%
  summarise(
    across(cut_val:pu, ~100*sd(.x, na.rm=T)),
    .groups='drop'
  ) %>%
  mutate(metric = 'sd')

# means
res_mean <- res_df %>%
  filter(trial == 0) %>%
  mutate(metric = 'mean') %>%
  select(colnames(res_sd))

# combined
res_df2 <- res_mean %>%
  bind_rows(res_sd) %>%
  pivot_wider(
    id_cols=group:sens_val,
    names_from=metric,
    values_from=cut_val:pu
  )

############
# Figure 1 #
############

fig1a <- res_df %>%
  filter(group == 'All participants') %>%
  ggplot(aes(x = 100*sens_val, y = 100*pu, color = predictor)) +
  geom_smooth(method = 'loess', aes(group=predictor), size=1.25, se=F) +
  facet_wrap(.~'Univariate plasma biomarkers') +
  ylim(c(0, 100)) +
  ylab('Tau PET scans saved (%)') +
  xlab('Biomarker sensitivity (%)') +
  theme_clean() +
  theme(legend.position = c(0.75,0.85), legend.title=element_blank(),
        legend.margin = margin(0., 0.2, 0.2, 0.2, "cm"),
        legend.background=element_rect(color='gray',
                                       size=0.5, linetype="solid")) +
  guides(color = guide_legend(nrow = 2))
fig1a <- ggpubr::set_palette(fig1a, 'jama')

fig1b <- res_df %>%
  filter(group == 'All participants') %>%
  filter(sens_val %in% c(0.90,0.95,0.99)) %>%
  ggplot(aes(x = factor(predictor, levels=rev(c('pTau217',
                                                'GFAP',
                                                'pTau181',
                                                'pTau231',
                                                'Ab42/Ab40',
                                                'NfL'))),
             y = 100*pu,
             color=predictor)) +
  geom_jitter(alpha=0.15) +
  stat_summary(fun = mean, geom = "crossbar",                                                                                                                                                                                                                                        size=0.5, width=0.7) +
  stat_summary(fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               size = 1, width=0.25,
               geom = "errorbar") +
  facet_wrap(.~paste0('Sensitivity: ', 100*sens_val, '%')) +
  theme_clean() +
  theme(legend.position='none', axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 335, vjust = 0., hjust=0.5)) +

  ylim(c(0,100)) +
  ylab('Tau PET scans saved (%)')
fig1b <- ggpubr::set_palette(fig1b, 'jama')


fig1 <- ggpubr::ggarrange(
  fig1a, fig1b,
  nrow=2,
  widths=c(0.6, 0.4), heights = c(0.6,0.4)
)

############
# Figure 2 #
############

fig3 <- res_df %>%
  filter(group == 'All participants') %>%
  filter(sens_val %in% c(0.90,0.95,0.975)) %>%
  ggplot(aes(y = fct_reorder(predictor, 100*ppv),
             x = 100*ppv,
             color=predictor)) +
  geom_vline(xintercept=38.4, linetype='dashed', color='gray', size=1) +
  stat_slab(aes(fill=predictor,
                     fill_ramp = stat(cut_cdf_qi(cdf, .width = c(0.02, 0.8, 0.95, 1)))),
            height = 1, color = "white", slab_size = 0.5) +
  scale_fill_viridis_d(option = "rocket", guide = "none", end = 0.9) +
  scale_fill_ramp_discrete(range = c(1, 0.2), guide = "none") +
  facet_wrap(.~paste0('Sensitivity: ', 100*sens_val, '%')) +
  theme_clean() +
  theme(legend.position='none', axis.title.y = element_blank()) +
  xlim(c(0,100)) +
  xlab('Positive predictive value (%)')

fig2 <- ggpubr::set_palette(fig2, 'jama')

# figure 3: pu at different sensitivity thresholds in subgroups relative to all
res_df0 <- res_df2 %>%
  filter(group == 'All participants') %>%
  select(-group)
colnames(res_df0)[5:ncol(res_df0)] <- paste0(colnames(res_df0)[5:ncol(res_df0)], '_ALL')
res_df3 <- res_df2 %>%
  filter(group != 'All participants') %>%
  left_join(res_df0, by=c('outcome','stat','predictor','sens_val')) %>%
  mutate(
    pu_mean_DIFF = pu_mean - pu_mean_ALL
  ) %>%
  filter(sens_val %in% c(0.9,0.95,0.99)) %>%
  mutate(
    pu_mean_DIFF = 100* pu_mean_DIFF,
    pu_mean_LO = pu_mean_DIFF - 1.96 * pu_sd,
    pu_mean_HI = pu_mean_DIFF + 1.96 * pu_sd,
    predictor_fct = factor(predictor,
                           levels=rev(c('pTau217',
                                    'GFAP',
                                    'pTau181',
                                    'pTau231',
                                    'Ab42/Ab40',
                                    'NfL'))),
    predictor = as.numeric(factor(predictor, levels=rev(c('pTau217',
                                                         'GFAP',
                                                         'pTau181',
                                                         'pTau231',
                                                         'Ab42/Ab40',
                                                         'NfL')))),
    group = factor(group, levels=c('SCD-MCI','All-cause dementia'))
  )

fig3 <- res_df3 %>%
  filter(
    sens_val %in% c(0.9, 0.95)
  ) %>%
  mutate(
    sens_val = paste0('Sensitivity: ', 100*sens_val, '%')
  ) %>%
  ggplot(aes(y = pu_mean_DIFF, x = predictor)) +
  geom_hline(yintercept=0, color='gray', linetype='dashed') +
  ggchicklet:::geom_rrect(aes(ymin = pu_mean_LO, ymax = pu_mean_HI,
                xmin=predictor-0.3, xmax=predictor+0.3, fill=group),
                r = unit(0.5, 'npc'), alpha = 0.5,
                position=position_dodge(width=0.9)
            ) +
  geom_point(aes(group=group, color=group, fill=group ),
             size=5, position=position_dodge(width=0.9)) +
  facet_nested_wrap(vars(group, sens_val), nrow = 2, strip = nested_settings) +
  theme_clean() +
  theme(axis.title.y = element_blank(), legend.position='none')+
  ylab('Difference in tau PET scans saved (%)') +
  scale_x_continuous(breaks=1:6, labels=levels(res_df3$predictor_fct)) +
  coord_flip()

fig3 <- ggpubr::set_palette(fig3, 'jama')


#ggsave('manuscript/Figure1.pdf', fig1, width = 12, height=9)
#ggsave('manuscript/Figure2.pdf', fig2, width = 12, height=4)
#ggsave('manuscript/Supp_FigureX.pdf', fig3, width = 8, height=6)




