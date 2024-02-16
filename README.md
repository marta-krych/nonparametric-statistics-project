
# Investigating Speech Patterns in Parkinson's Disease (PD) and  Rapid Eye Movement Sleep Behavior Disorder (RBD)

This project is inspired by a seminal study [1] that compared speech recordings from individuals with Rapid Eye Movement Sleep Behavior Disorder (RBD), newly diagnosed PD patients, and healthy controls. 
This study demonstrated the efficacy of automated vocal analysis in detecting subliminal parkinsonian speech deficits in RBD patients, a group at high risk for developing PD. 
The findings underscore the potential of automated methods in identifying early indicators of neurodegeneration, paving the way for timely interventions.
The implications of automated vocal analysis for the screening and diagnosis of neurodegenerative disorders are profound. By facilitating early detection of PD and other synucleinopathies, this approach could significantly contribute to the development of targeted therapies and improve the quality of life for affected individuals. Our project aims to build upon these findings by applying nonparametric statistical methods to enhance the understanding and detection of neurodegenerative patterns in speech, thereby contributing to the evolving landscape of PD diagnosis and management.

### Dataset

Hlavnika,J., Tykalov,T., Onka,K., Rika,E., Rusz,J., and J.,J.. (2017). Early biomarkers of Parkinson’s disease based on natural connected speech. UCI Machine Learning Repository. https://doi.org/10.24432/C5W02Q.


### Analysis

We applied nonparametric methods for feature selection (ANOVA/MANOVA, clustering and permutation tests) to the speech-related variables in order to build a logistic regression model. 
The aim of the model is to assign to each RBD patients a probability of beign classified as Parkinson ill (to be interpreted as the risk of developing the disease).


### Authors

Ange Dakouri,
Simone Giacomello,
Marta Krychkovska,
Antonio Napolitano

### Institution
Nonparametric Statistics project for the academic year 2023/2024 @PoliMi.

### Acknowledgments

[1] Hlavnička, J., Čmejla, R., Tykalová, T. et al. Automated analysis of connected speech reveals early biomarkers of Parkinson’s disease in patients with rapid eye movement sleep behaviour disorder. Sci Rep 7, 12 (2017). https://doi.org/10.1038/s41598-017-00047-5

