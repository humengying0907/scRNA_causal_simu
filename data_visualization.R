
library(dplyr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(fitdistrplus)

######### prepare input for FR-Perturb (use co-culture for now) ##############
meta = read.delim('SCP1064/meta.csv',sep = ',')
meta_coculture = meta[meta$condition == 'Co-culture',]

all_targets = meta_coculture %>%
  filter(!is.na(targets)) %>%                # Remove NA values
  separate_rows(targets, sep = ",") %>%      # Split and expand rows on comma
  distinct(targets)  
all_targets = all_targets$targets
all_targets = all_targets[all_targets!='']

target_matrix <- matrix(0, ncol = length(meta_coculture$NAME), nrow = length(all_targets)+1, 
                        dimnames = list(c(all_targets,'unperturbed'),c(meta_coculture$NAME)))


for(i in 1:ncol(target_matrix)){
  targets = strsplit(as.character(meta_coculture$targets[i]), ",")[[1]]
  if(length(targets)>0){
    for (target in targets) {
      target_matrix[target,i] <- 1
    }
  }else{
    target_matrix['unperturbed',i] = 1
  }
}

write.table(target_matrix,file = 'FR_Perturb/coculture_perturbation_matrix.txt',sep = '\t',quote = FALSE)


################### check FR results ##############
# load effect size table
coculture_lfc = read.delim('FR_Perturb/cocult_LFCs.txt',sep = '\t',row.names = 1)

sgRNA_candidate = c('STAT1','JAK1','JAK2','IFNGR1','IFNGR2')
downstream_genes = c('PSMA3', 'PSMA1', 'PSMB1', 'B2M', 'PSMA7', 'PSMB8', # Antigen processing and presentation of exogenous peptide antigen via MHC class I
                     'CENPE', 'HLA-DMA', 'RACGAP1', 'KIF4A', 'HLA-DPB1', 'HLA-DRA', 'KIF23', 'KIF2C', 'HLA-DRB1', 'HLA-DPA1' #Antigen processing and presentation of exogenous peptide antigen via MHC class II
                     )

mat = coculture_lfc[downstream_genes,sgRNA_candidate]
Heatmap(mat,col = colorRamp2( seq(-5,0,0.25), colors = hcl.colors(21, "Blues")),rect_gp = gpar(col = "#a9a9a9", lwd = 0.5),name = 'lfc')

################### explore NB paramters #################
# explore scExpr data and check what is the parameter of NB
ad <- anndata::read_h5ad('SCP1064/cocult/raw_unperturbed_adata.h5ad')
meta = ad$obs
ad$X[1:5,1:5]

# IFN-gamma reponse pathway related genes
hist(ad$X[,'STAT1'])
fit <- fitdistrplus::fitdist(ad$X[,'STAT1'], "nbinom")
hist(rnbinom(1000, size = 1.54, mu = 7.48))
#p = 1.54/(1.54+7.48)
#hist(rnbinom(1000, size = 1.54, p = p)) # size is equivalent to disperison parameter in scipy.stats.nbinom.rvs
print(fit$estimate)
# size       mu 
# 1.546387 7.491518 



hist(ad$X[,'IRF1'])
fit <- MASS::fitdist(ad$X[,'IRF1'], "nbinom")
hist(rnbinom(1000, size = 1.449556, mu = 3.629281))
print(fit$estimate)
# size       mu 
# 1.449556 3.629281 

############### explore Gaussian parameters #################
ad <- anndata::read_h5ad('SCP1064/cocult/gene_filtered_adata.h5ad')
meta = ad$obs
ad_unperturbed <- ad[ad$obs$MOI == 0,]
# anndata::write_h5ad(ad_unperturbed,filename = 'SCP1064/cocult/log1p_unperturbed_adata.h5ad')

fit = MASS::fitdistr(ad_unperturbed$X[,'STAT1'], "normal") # mean 3.722293415 ; sd 1.197740901
fit = MASS::fitdistr(ad_unperturbed$X[,'STAT1'][ad_unperturbed$X[,'STAT1']!=0], "normal") # mean 3.987191528 ; sd 0.693163974


fit = MASS::fitdistr(ad_unperturbed$X[,'IRF1'], "normal") # mean 2.84831820 ; sd 1.40570746
fit = fitdistrplus::fitdist(ad_unperturbed$X[,'IRF1'][ad_unperturbed$X[,'IRF1']!=0], "norm") # mean 3.388740716 ; sd 0.720817872



################## evaluate simulation results ################
Gaussian_raw = read.delim('project_dir/Gaussian_raw_performance.txt',sep = '\t')
NB_raw = read.delim('project_dir/NB_raw_dispersion_2_performance.txt',sep = '\t')
# NB_log = read.delim('project_dir/NB_log_dispersion_2_performance.txt',sep = '\t')



# SHD in Z only
Z_performance = rbind(Gaussian_raw[Gaussian_raw$data_type == 'Oracle',] %>% mutate(group = 'Gaussian'),
                      NB_raw[NB_raw$data_type == 'Oracle',] %>% mutate(group = 'NB'),
                      NB_log[NB_log$data_type == 'Oracle',] %>% mutate(group = 'NB_log'))

Z_performance_summary = Z_performance %>% group_by(node, degree,test_method,group) %>% 
  summarise(avg_SHD = mean(SHD, na.rm = T),
            max_SHD = max(SHD, na.rm = T),
            min_SHD = min(SHD, na.rm = T))
# Z_performance_summary = Z_performance_summary[Z_performance_summary$node %in% c(10,20),]
Z_performance_summary$node =paste0('node = ',Z_performance_summary$node)
Z_performance_summary$node = factor(Z_performance_summary$node,levels =unique(Z_performance_summary$node) )
#Z_performance_summary = Z_performance_summary[Z_performance_summary$group!='Gaussian',]
Z_performance_summary = Z_performance_summary[Z_performance_summary$group!='NB_log',]
Z_performance_summary$group[Z_performance_summary$group == 'NB'] = 'Negative Binomial'
#Z_performance_summary$group[Z_performance_summary$group == 'NB_log'] = 'Log transformed Negative Binomial'
Z_performance_summary = Z_performance_summary[Z_performance_summary$test_method == 'zerodel_fisherz',]

pd <- position_dodge(width = 0.1)

ggplot(Z_performance_summary,aes(x = degree, 
                                     y = avg_SHD, group = group)) +
  geom_line(position = pd,aes(color=group)) +
  geom_errorbar(aes(ymin = min_SHD, ymax = max_SHD, color=group),width = .1, position = pd, linetype = 1)+
  geom_point(position = pd,aes(color=group )) +
  # scale_color_manual(values = c('#ffab4e','#00BFC4'))+
  guides(linetype = guide_legend("Transmission"))+
  facet_grid(test_method~node)+
  ylab('SHD')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        strip.background = element_rect(fill = 'white'),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9))+
  guides(color=guide_legend(title="True expression Z"))


ggplot(Z_performance_summary,aes(x = degree, 
                                 y = avg_SHD, group = test_method)) +
  geom_line(position = pd,aes(color=test_method)) +
  geom_errorbar(aes(ymin = min_SHD, ymax = max_SHD, color=test_method),width = .1, position = pd, linetype = 1)+
  geom_point(position = pd,aes(color=test_method )) +
  # scale_color_manual(values = c('#ffab4e','#00BFC4'))+
  guides(linetype = guide_legend("Transmission"))+
  facet_grid(group~node)+
  ylab('SHD')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        strip.background = element_rect(fill = 'white'),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9))+
  guides(color=guide_legend(title="PC method"))

# NB vs Gaussain vs Gaussain = 0.6 dropout
Z_performance = rbind(Gaussian_raw[Gaussian_raw$data_type == 'Oracle',] %>% mutate(group = 'Gaussian'),
                      NB_raw[NB_raw$data_type == 'Oracle',] %>% mutate(group = 'Negative Binomial'),
                      Gaussian_raw[Gaussian_raw$drop_rate_levels == '0.6',] %>% mutate(group = 'Gaussian with 0.6 dropout'))

Z_performance_summary2 = Z_performance %>% group_by(degree,node,test_method,group) %>%
  summarise(avg_SHD = mean(SHD, na.rm = T),
            max_SHD = max(SHD, na.rm = T),
            min_SHD = min(SHD, na.rm = T))
Z_performance_summary2 = Z_performance_summary2[Z_performance_summary2$test_method == 'zerodel_fisherz',]
Z_performance_summary2$node =paste0('node = ',Z_performance_summary2$node)
Z_performance_summary2$node = factor(Z_performance_summary2$node,levels =levels(Z_performance_summary$node) )

ggplot(Z_performance_summary2,aes(x = degree, 
                                 y = avg_SHD, group = group)) +
  geom_line(position = pd,aes(color=group)) +
  #geom_errorbar(aes(ymin = min_SHD, ymax = max_SHD, color=group),width = .1, position = pd, linetype = 1)+
  geom_point(position = pd,aes(color=group )) +
  scale_color_manual(values = c("#F8766D" ,'#ffab4e',"#00BFC4"))+
  guides(linetype = guide_legend("Transmission"))+
  facet_grid(test_method~node)+
  ylab('SHD')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        strip.background = element_rect(fill = 'white'),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9))+
  guides(color=guide_legend(title="simulated expression"))


# effect of dropout rate 
add_droput_cut = function(df){
  df$drop_rate_levels = cut(df$dropout_rate, breaks = c(-Inf, 0.001, 0.3, 0.5,Inf),
                            right = FALSE,labels = c("0", "0.2", "0.4", "0.6"))
  return(df)
}

Gaussian_raw = add_droput_cut(Gaussian_raw)
NB_raw = add_droput_cut(NB_raw)
NB_log = add_droput_cut(NB_log)

ggplot(NB_raw,aes(drop_rate_levels,SHD))+
  geom_boxplot(aes(color = test_method))+
  facet_grid(degree~exo_noise_mu)+
  theme_bw()

Gaussian_raw_dropout_summary = Gaussian_raw[Gaussian_raw$node == 20,] %>% 
  group_by(degree,test_method,drop_rate_levels) %>%
  summarise(avg_SHD = mean(SHD, na.rm = T),
            max_SHD = max(SHD, na.rm = T),
            min_SHD = min(SHD, na.rm = T))

Gaussian_raw_dropout_summary$degree = as.factor(Gaussian_raw_dropout_summary$degree)
ggplot(Gaussian_raw_dropout_summary, aes(x = drop_rate_levels, y = avg_SHD,group = degree)) +
  geom_line(aes(color = degree), position = pd) +
  geom_point(aes(color = degree), position = pd) +
  geom_errorbar(aes(ymin = min_SHD, ymax = max_SHD, color = degree), width = 0.1, position = pd) +
  guides(linetype = guide_legend("Transmission"))+
  facet_wrap(~test_method, nrow = 1)+
  ylab('SHD')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        strip.background = element_rect(fill = 'white'),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9))


# effect of exo_noise_mu
mu_summary = rbind(Gaussian_raw[Gaussian_raw$drop_rate_levels == '0.6',] %>% mutate(group = 'Gaussian'),
                   NB_raw[NB_raw$data_type == 'Oracle',] %>% mutate(group = 'Negative Binomial'))
mu_summary = mu_summary[mu_summary$degree != 5,]
mu_summary$degree = as.factor(mu_summary$degree)
mu_summary$exo_noise_mu = as.factor(mu_summary$exo_noise_mu)
mu_summary = mu_summary[mu_summary$test_method == 'zerodel_fisherz',]
# mu_summary$node =paste0('node = ',mu_summary$node)
# mu_summary$node = factor(mu_summary$node,levels = levels(Z_performance_summary$node))


ggplot(mu_summary,aes(degree,SHD))+
  geom_boxplot(aes(color = exo_noise_mu))+
  facet_grid(node~group,scales = 'free_x')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        strip.background = element_rect(fill = 'white'),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9))



# check simulated expression distribution
library(reticulate)
np <- import("numpy")

Z_Gaussian_10_node_3_degree = np$load("project_dir/simulated_data/Z_Gaussian_exo_noise_mu_3.0_rep_0.npy")
Z_NB_10_node_3_degree = np$load("project_dir/simulated_data/Z_NB_exo_noise_mu_3.0_rep_0.npy")

par(mfrow = c(1, 2))
hist(Z_Gaussian_10_node_3_degree[,3],xlab = 'simulated expression value',main = 'Z-Gaussian')
hist(Z_NB_10_node_3_degree[,3],breaks = 20,xlab = 'simulated expression value',main = 'Z-Negative binomial')


