# Data and Code for the paper "Hyper-cores promote localization and efficient seeding in higher-order processes"
This repository contains the data and code associated to the paper "Hyper-cores promote localization and efficient seeding in higher-order processes" by Marco Mancastroppa, Iacopo Iacopini, Giovanni Petri and Alain Barrat, [arXiv:2301.04235](https://arxiv.org/abs/2301.04235) (2023)
# Data
The data that support the findings of this study are publicly available:
* **SocioPattern** data sets (InVS15, LH10, SFHH, LyonSchool, Thiers13) by the [SocioPatterns project](http://www.sociopatterns.org/). Data source [here](http://www.sociopatterns.org/datasets/);
* **Utah’s schools** data sets (Mid1, Elem1) by the Contacts among Utah's School-age Population (CUSP), presented in [Toth et al. J. R. Soc. Interface 12: 20150279 (2015)](https://royalsocietypublishing.org/doi/10.1098/rsif.2015.0279). Data source [here](https://royalsocietypublishing.org/doi/suppl/10.1098/rsif.2015.0279);
* **Email-EU** data set, presented in [A. Paranjape et al. Proceedings of the Tenth ACM International Conference on Web Search and Data Mining, p. 601–610 (2017)](https://dl.acm.org/doi/10.1145/3018661.3018731). Data source [here](https://www.cs.cornell.edu/~arb/data/)
* **Email-Enron** data set, presented in [A. R. Benson et al., PNAS 115, E11221 (2018)](https://www.pnas.org/doi/10.1073/pnas.1800683115). Data source [here](https://www.cs.cornell.edu/~arb/data/)
* **Political interactions** data sets (congress-bills, senate-bills, house-committees, senate-committees), presented in [P. S. Chodrow et al., Science Advances 7, eabh1303 (2021)](https://www.science.org/doi/10.1126/sciadv.abh1303), C. Stewart III et al., Congressional committee assignments, 103rd to 114th congresses (1993–2017), [J. H. Fowler, Social Networks 28, 454 (2006)](https://doi.org/10.1016/j.socnet.2005.11.003) and [J. H. Fowler, Political Analysis 14, 456–487 (2006)](https://doi.org/10.1093/pan/mpl002). Data source [here](https://www.cs.cornell.edu/~arb/data/)
* **Online interactions** data sets (music-review, algebra-questions, geometry-questions), presented in [J. Ni et al., Proceedings of the 2019 EMNLP-IJCNLP, pp. 188–197 (2019)](https://aclanthology.org/D19-1018/) and [I. Amburg et al., Proceedings of the 2022 SIAM International Conference on Data Mining (SDM), pp. 145–153 (2022)](https://epubs.siam.org/doi/10.1137/1.9781611977172.17). Data source [here](https://www.cs.cornell.edu/~arb/data/)
* **Ecological** data sets (M_PL_015_ins, M_PL_015_pl, M_PL_062_ins, M_PL_062_pl) by the [Web of life: ecological network database](https://www.web-of-life.es ). Data source [here](https://www.web-of-life.es).

Both the original and processed data are collected in the `Data` folder, with the preprocessing code to obtain the empirical static hypergraphs used in the study.
# Code

The `Hyper-core_decomposition` folder contains the code to obtain the (k,m)-hyper-core decomposition of a static hypergraph, and also the k-core and s-core decomposition of the associated projected graph. 

The `Hypergraph_randomization` folder contains the code to obtain a randomized realization of the original static hypergraph, through hyperedge reshuffling. The code is an adaptation of the reshuffling procedure proposed in [N. W. Landry et al., Chaos: An Interdisciplinary Journal of Nonlinear Science 32, 053113 (2022)](https://doi.org/10.1063/5.0086905) (the original procedure can be found at [https://github.com/nwlandry/hypergraph-assortativity](https://github.com/nwlandry/hypergraph-assortativity)).  

The `Nonlinear_higher-order_contagion` folder contains the code to simulate the SIS and SIR higher-order nonlinear contagion processes on static hypergraphs.  

The `Threshold_higher-order_contagion` folder contains the code to simulate the SIS and SIR higher-order threshold contagion processes on static hypergraphs.  

The code to simulate the naming-game process on hypergraphs with committed minority is available at [https://github.com/iaciac/higher-order-NG](https://github.com/iaciac/higher-order-NG)

The code uses the CompleX Group Interactions (XGI) library in Python [https://xgi.readthedocs.io](https://xgi.readthedocs.io/)  
XGI repository: [https://github.com/xgi-org/xgi](https://github.com/xgi-org/xgi)  
Landry, N. W., Lucas, M., Iacopini, I., Petri, G., Schwarze, A., Patania, A., & Torres, L. (2023). XGI: A Python package for higher-order interaction networks. Journal of Open Source Software, 8(85), 5162. [https://doi.org/10.21105/joss.05162](https://doi.org/10.21105/joss.05162)
