# ZICRM-ASCA+
Zero-inflated Counts RM-ASCA+
_Current advancements in biological data generation require comparative analysis methods which account for the design of the experiment and the specific distributions underlying the data. In microbiome analysis, the presence of marker gene sequences (zOTUs) are used to quantify microbial compositions. These, in turn, can be used for assessment of treatment-specific responses. As experimental designs can be complex, with longitudinal analysis of microbial changes in multiple, unbalanced treatment groups, a method is required for accurate determination of the effects of these different factors. For data with a Gaussian distribution, repeated measures-ANOVA Simultaneous Component Analysis+ (RM-ASCA+) provides such a method, combining a generalized linear mixed model (GLMM) decomposition and principal component analysis (PCA) for visualization of latent structure in design factor effects. However, RM-ASCA+ is not suited for overdispersed, zero-inflated count data, which is generally modelled with a negative binomial (NB) distribution. Here we present a zero-inflated counts RM-ASCA+ (ZICRM-ASCA+) extension for comparative microbiome analysis, which models zOTU counts with a zero-inflated NB mixed model. A simulation study is used to show the performance of ZICRM-ASCA+ in capturing latent structure in experimental design effects for zero-inflated, overdispersed data. Furthermore, the method is applied to the gut microbiome data set of the  CALM intervention study, which reveals shared effects between zOTUs with a high and low level of zero-inflation, as well as an unexpectedly high between-group variation at baseline. ZICRM-ASCA+ builds upon a framework for multivariate data analysis to provide an exploratory analysis of latent structure which, alone or combined with other quantitative data, could be used to reveal new connections between biological features._
- Abstract from Master's Thesis 

This GitHub page is dedicated to the simulation study for ZICRM-ASCA+. The methodology is based on RM-ASCA+ (Madssen et al. 2020) and glmmTMB (Brooks et al. 2017) among other studies. The simulation can simply be run by downloading the entire project and running the RMarkdown file in Rstudio.

If you encounter any problems, feel free to let me know. 

Auke



Brooks,  M.  E.,  Kristensen,  K.,  van  Benthem,K.  J.,  Magnusson,  A.,  Berg,  C.  W.,Nielsen, A., Skaug, H. J., Maechler, M.,& Bolker, B. M. (2017). glmmTMB bal-ances speed and flexibility among pack-ages for zero-inflated generalized linearmixed modeling.The  R  Journal,9(2). 

Madssen, T. S., Giske deg ard, G. F., Smilde, A. K., & Westerhuis, J. A. (2020). Repeated measures ASCA+ for analysis of longitudinal intervention studies with multivariate outcome data. medRxiv. https://doi.org/10.1101/2020.12.03.20243097

