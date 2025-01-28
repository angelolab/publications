# SOTIP  v1.2

## SOTIP is a Versatile Method for Microenvironment Modelling with Spatial Omics Data

### Developer: Zhiyuan Yuan (zhiyuan AT fudan DOT edu DOT cn) Yisi Li (li-ys16 AT mails DOT tsinghua DOT edu DOT cn)

The rapid development of spatial omics techniques generates datasets with diverse scales and modalities. However, most existing methods focus on modeling dynamics of single cells and ignore microenvironments (MEs) which bridge single cells to tissues. Here we present SOTIP, a scalable framework incorporating MEs and their relationships into a unified graph. Based on this graph, three tasks can be performed, namely, spatial heterogeneity (SHN) quantification, spatial domain (SDM) identification and differential microenvironment (DME) analysis. We validate SOTIP’s high performance of accuracy, robustness, interpretation, and scalability on various datasets by comparing with state-of-art methods. In mammalian cerebral cortex, we reveal a striking gradient SHN pattern with strong correlations with the cortical depth. In human triple-negative breast cancer (TNBC), we identify previously unreported molecular polarizations around SOTIP-detected tumor-immune boundaries. Most importantly, we discover MEs which specifically occur in different TNBC subtypes with certain compositional and spatial properties, could be powerful clinical indicators. Overall, by modeling biologically explainable microenvironments, SOTIP outperforms state-of-art methods and provides new perspectives for data interpretation, which facilitates further understanding of spatial information on a variety of biological issues.

![SOTIP workflow](SOTIP_analysis/images/workflow.png)


## Usage

The [**SOTIP**](https://github.com/TencentAILabHealthcare/SOTIP) package is a graph based analytical framework for various spatial omics data. There are three modules for SOTIP:

- Spatial heterogeneity (SHN) quantification. [**4i**](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/tutorial/4i_HeLa.ipynb),  [**Visium**](https://github.com/yuanzhiyuan/SOTIP/tree/master/SOTIP_analysis/tutorial/Visium_Zebrafish.ipynb)
- Spatial domain (SDM) identification. [**MIBI**](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/tutorial/MIBI_TNBC.ipynb),  [**scMEP**](https://github.com/yuanzhiyuan/SOTIP/tree/master/SOTIP_analysis/tutorial/scMEP_CLCC.ipynb),  [**Visium**](https://github.com/yuanzhiyuan/SOTIP/tree/master/SOTIP_analysis/tutorial/Visium_Cortex.ipynb),  [**osmFISH**](https://github.com/yuanzhiyuan/SOTIP/tree/master/SOTIP_analysis/tutorial/osmFISH_cortex.ipynb)
- Differential microenvironment (DME) analysis. [**MIBI**](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/DMA_TNBC.ipynb)
- SOTIP supports various spatial omics protocols ranging from spatial transcriptomics (10x Visium, osmFISH, seqFISH+), spatial proteomics (4i, MIBI-TOF, scMEP), and spatial metabolomics (TOF-SIMS).

#### Example of SDM identification on 10x Visium

<p align="center">
<img src="SOTIP_analysis/images/merge_Visium.gif"/>
<br>
Herarchical merging for 10x Visium data
</p>

#### Example of SDM identification on MIBI-TOF

<p align="center">
<img src="SOTIP_analysis/images/merge_MIBI.gif"/>
<br>
Herarchical merging for MIBI-TOF data
</p>

## Tutorial
Please install Jupyter in order to open these notebooks.

For the step-by-step tutorial, please refer to: 
<br>
https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/tutorial/
<br>

For the reproduction of paper's results, please refer to:
<br>
https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/
<br>

Please download demo datasets from following doi: 
<br>
10.6084/m9.figshare.18516128.
<br>

## How to install?
- git clone this repository
- python setup.py install
- conda install pyemd

## SOTIP has been tested on

- System: CentOS
- Python: 3.8.0
- Python packages: numpy==1.21.2 pandas==1.3.4 scipy==1.7.1 matplotlib==3.4.3 seaborn==0.11.2 scanpy==1.8.2  squidpy==1.1.2 palettable==3.3.0 scikit-learn==1.0.1 networkx==2.6.3 shapely==1.8.0 pyemd==0.5.1

## Reproduce

|  Figure |Description   |  Notebook |
| :------------: | :------------: | :------------: |
|  Fig 2a-b | Simulation data1  |  [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/simulation/simulation1.ipynb) |
| Fig 2c-d  | Simulation data2  | [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/simulation/simulation2.ipynb)  |
|  Fig 2e-f | Simulation data3  | [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/simulation/simulation3.ipynb)  |
|Fig 3a-f   |SOTIP-SNH on 4i data   | [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/4i/)  |
|  Fig 3g-i |  SOTIP-SHN on Zebrafish Visium data  |  [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/Visium_Zebrafish/) |
| Fig 3k,n  |SOTIP-SHN on Mouse cortex osmFISH data   | [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/osmFISH_Cortex/)  |
| Fig 3m,o  | SOTIP-SHN on Mouse cortex seqFISH+ data  |  [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/seqFISH_Cortex/) |
| Fig 4b-c  |SOTIP-SDM on SpatialLIBD Visium data   |  [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/Visium_Cortex/SDM_Visium_cortex.ipynb) |
| Fig 4d  | SOTIP-SDM on SpatialLIBD Visium benchmark  |  [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/Visium_Cortex/SpatialLIBD_stats.ipynb) |
|  Fig 4g-h | SOTIP-SDM on Mouse cortex osmFISH benchmark  | [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/osmFISH_Cortex/SDM_osmFISH.ipynb)  |
| Fig 5e (right)  | SOTIP-SDM on Human CRC scMEP data  |  [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/scMEP_CLCC/sotip_result/) |
|  Fig 5e (left) | SOTIP-SDM on Human CRC scMEP data  | [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/scMEP_CLCC/spagcn_result/)  |
| Fig 5f   | SOTIP-SDM on Human CRC scMEP data  |  [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/scMEP_CLCC/) |
|  Fig 5h (top) |  SOTIP-SDM on Human TNBC MIBI data | [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/TNBC_p4_cur_20211212.ipynb)  |
| Fig 5h (bottom)  |  SOTIP-SDM on Human TNBC MIBI data|  [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/TNBC_p9_cur_20211212.ipynb) |
|Fig 5i (top)   | SOTIP-SDM on Human TNBC MIBI data  |  [link1](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/compare_polar_rst_p4.ipynb)  [link2](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/sotip_MIBI_polarization_p4_rstest.ipynb) [link3](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/spagcn_MIBI_polarization_p4_rstest.ipynb)|
| Fig 5i (bottom)  | SOTIP-SDM on Human TNBC MIBI data  | [link1](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/compare_polar_rst_p9.ipynb)  [link2](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/sotip_MIBI_polarization_p9_rstest.ipynb) [link3](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/spagcn_MIBI_polarization_p9_rstest.ipynb)|
| Fig 6  | SOTIP-DME on SIMS liver data  | [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/SIMS_liver/)  |
| Fig 7a-j  | SOTIP-DME on MIBI TNBC data  |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/DMA_TNBC.ipynb)|
| Fig 7k-m  | Validation of SOTIP-DME on MIBI TNBC data  |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/validate_DMA_save.ipynb)|
|  Fig 7n |   Validation of TNBC and survival analysis | [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/validate_DMA_load.ipynb)|
| SuppFig 1a-c  | Simulation data 1  |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/simulation/simulation1_stats.ipynb)|
|  SuppFig 1d-g | Simulation data 2  |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/simulation/simulation2_stats.ipynb)|
|  SuppFig 1h-k | Simulation data 3  |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/simulation/simulation3_stats.ipynb)|
| SuppFig 2a-b  |  Mouse cortex osmFISH data |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/osmFISH_Cortex/)|
|  SuppFig 2c-d |   Mouse cortex seqFISH+ data|   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/seqFISH_Cortex/)|
| SuppFig 2e-f  | Zebrafish Visium data   |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/Visium_Zebrafish/SHN_Visium_zebroFISH_A.ipynb)|
|  SuppFig 4g |  Human TNBC MIBI data |   [link1](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/compare_polar_rst_p9.ipynb)  [link2](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/sotip_MIBI_polarization_p9_rstest.ipynb) [link3](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/spagcn_MIBI_polarization_p9_rstest.ipynb)|
| SuppFig 5  | SOTIP-DME on SIMS liver data  |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/SIMS_liver/)|
| SuppFig 6  | SOTIP-DME on MIBI TNBC data  |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/DMA_TNBC.ipynb)|
| SuppFig 7  | Validation of TNBC and survival analysis  |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/MIBI_TNBC/validate_DMA_load.ipynb)|
| SuppFig 8  | Simulation data 4  |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/Simulation/simulation4.ipynb)|
|  SuppFig 11 |   Simulation data 5|   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/Simulation/simulation5.ipynb)|
| SuppFig 12  |  Simulation data 6 |   [link](https://github.com/TencentAILabHealthcare/SOTIP/tree/master/SOTIP_analysis/Simulation/simulation6.ipynb)|







## Disclaimer

This tool is for research purpose and not approved for clinical use.

This is not an official Tencent product.

## Copyright

This tool is developed in Tencent AI Lab.

The copyright holder for this project is Tencent AI Lab.

All rights reserved.


## References

Please consider citing the following reference:

- Yuan, Z., Li, Y., Shi, M. et al. SOTIP is a versatile method for microenvironment modeling with spatial omics data. Nat Commun 13, 7330 (2022). https://doi.org/10.1038/s41467-022-34867-5

<br>

