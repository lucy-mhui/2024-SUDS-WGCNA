# Abstract

Major depressive disorder (MDD) is the most prevalent psychiatric condition marked by persistent sadness and cognitive impairments. The monoamine hypothesis, which attributes depressive symptomatology to dysfunctional neurotransmitter systems, has formed the basis for contemporary treatments. However, only 50% of the population respond to the treatments, implicating alternative pathophysiological underpinnings. Recent research indicates that neuroinflammatory processes play a significant role in its development. Importantly, there are notable sex differences in both the presentation and underlying neuroinflammatory mechanisms. However, many studies on this topic have conducted univariate analyses of inflammatory biomarkers, which fails to capture multidimensional interactions. A well-established network modeling technique called Weighted Gene Co-expression Network Analysis (WGCNA) may unveil gene co-expression modules in MDD. 

Using the Canadian Biomarker Integration Network in Depression (CAN-BIND-1) dataset, this project aimed to identify sex-dependent modules of co-expressed genes associated with inflammatory biomarkers in MDD patients and controls. Furthermore, this project investigated the correlation between inflammatory gene expression modules, clinical variables, and neuroimaging markers to elucidate the clinical relevance of immune dysregulation in MDD. Biomarkers from 211 MDD participants and 122 controls were included across three modalities: molecular, clinical, and neuroimaging. Blood samples were collected for the measurement of 29 chemokine/cytokine levels. WGCNA was employed for identifying disease-specific and sex-specific genetic clusters. Correlation between genetic modules and major depressive disorder was evaluated using clinical and neuroimaging measurements. 

About 76% of the genes were found to form five co-expression modules in the MDD cohort. Two out of the five modules were significantly correlated with the Montgomery-Asberg Depression Rating Scale (MADRS) (p<0.05). Of note, a network of 3 eigengenes (interleukin 17, fibroblast growth factor, and platelet-derived growth factor) exhibited the highest correlation to MADRS (r2=0.26). Interestingly, modules identified in the MDD cohort were negatively correlated with the female sex and positively correlated with the male sex, while the opposite was observed for the controls. Future results will include module correlations with diffusion weighted imaging markers across both cohorts. In summary, our results corroborate previous research on molecular-level immune dysfunction in MDD and offer new perspectives on the underlying mechanisms of the disorder.

<img width="404" alt="image" src="https://github.com/user-attachments/assets/3b185a98-beef-4671-a92d-f0d72daaa522">

# CAN-BIND-1 data structure

<img width="413" alt="image" src="https://github.com/user-attachments/assets/ef6a2ccb-4247-42d6-9de0-f2ac2f3c384a">

# Github code structure
```
├── 1 inflammatory markers             # Documentation files (alternatively `doc`)
│   ├── 1.1 descriptive statistics     # Table of contents
│   ├── 1.2 data input                 # Frequently asked questions
├── 2 depression                       # Documentation files (alternatively `doc`)
├── 3 dMRI                             # Documentation files (alternatively `doc`)
└── ...
```

# Methods
* 

# References and resources
1. Miller, A., Raison, C. The role of inflammation in depression: from evolutionary imperative to modern treatment target. Nat Rev Immunol 16, 22–34 (2016). https://doi.org/10.1038/nri.2015.5
2. Horvath WGCNA tutorial: https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=2&dl=0
3. CAN-BIND-1 description: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6606427/
