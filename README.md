# Hyper-cores promote localization and efficient seeding in higher-order processes
This repository contains the code associated to the paper "Hyper-cores promote localization and efficient seeding in higher-order processes" by Marco Mancastroppa, Iacopo Iacopini, Giovanni Petri and Alain Barrat, [arXiv:2301.04235](https://arxiv.org/abs/2301.04235) (2023)
# Data
The data that support the findings of this study are publicly available:
* **SocioPattern** data sets (InVS15, LH10, SFHH, LyonSchool, Thiers13) by the [SocioPatterns project](http://www.sociopatterns.org/). Data source [here](http://www.sociopatterns.org/datasets/);
* **Utah’s schools** data sets (Mid1, Elem1) by the Contacts among Utah's School-age Population (CUSP), presented in [Toth et al. J. R. Soc. Interface 12: 20150279 (2015)(https://royalsocietypublishing.org/doi/10.1098/rsif.2015.0279}. Data source [here](https://royalsocietypublishing.org/doi/suppl/10.1098/rsif.2015.0279);
* **Email-EU** data set. Data source [here](https://www.cs.cornell.edu/~arb/data/)
* **Email-Enron** data set. Data source [here](https://www.cs.cornell.edu/~arb/data/)
* **Political interactions** data sets. Data source [here](https://www.cs.cornell.edu/~arb/data/)
* **Online interactions** data sets. Data source [here](https://www.cs.cornell.edu/~arb/data/)
* **Ecological** data sets by the Web of life ecological. Data source [here](https://www.web-of-life.es).

Both the original and processed data are collected stored in the Data folder, with the preprocessing code to obtain the corresponding static hypergraphs.
# Code

The code used to simulate the Naming-Game (NG) process on hypergraphs with committed minority is available at [https://github.com/iaciac/higher-order-NG](https://github.com/iaciac/higher-order-NG)

The code use the CompleX Group Interactions (XGI) library in Python [https://xgi.readthedocs.io](https://xgi.readthedocs.io/)  
XGI repository: [https://github.com/xgi-org/xgi](https://github.com/xgi-org/xgi)  
Landry, N. W., Lucas, M., Iacopini, I., Petri, G., Schwarze, A., Patania, A., & Torres, L. (2023). XGI: A Python package for higher-order interaction networks. Journal of Open Source Software, 8(85), 5162. [https://doi.org/10.21105/joss.05162](https://doi.org/10.21105/joss.05162)
