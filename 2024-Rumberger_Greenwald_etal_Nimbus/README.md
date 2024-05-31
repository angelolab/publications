# Code for Rumberger, Greenwald et al., Improving cell phenotyping by factoring in spatial marker expression patterns with Nimbus

### User-friendly Jupyter notebooks for running Nimbus are available at [github.com/angelolab/Nimbus-Inference](https://github.com/angelolab/Nimbus-Inference). Please see that repository to use Nimbus for your own data.

The scripts here were used to generate the figures in the paper. Description of files:
1. ct_assignment.py: Assigns granular cell types from individual publications to the 10 lineages used in the paper.
2. figure_2.ipynb: Contains code for generating the overview of the Pan-M dataset.
    - Figure 2a-i
3. figure_3.ipynb: Contains code for generating the qualitative view on the predictions.
    - Figure 3a-e
4. figure_4.ipynb: Contains code for generating the quantitative view on the predictions.
    - Figure 4a-c and Supplement Figure 1h
5. figure_4_error_analysis_supp2.ipynb: Contains code for generating the error analysis.
    - Figure 4d and Supplement Figure 2 a-e
6. figure_5.ipynb: Contains code for generating cell phenotypes using Nimbus
    - Figure 5a-e
7. figures_supplement.ipynb: Contains additional supplementary figures
    - Supplement Figure 1a-g and Supplement Figure 2f

The notebooks expect the following folder structure to be present in the same directory as the notebooks:
```bash
.
├── data
│   ├── DCIS
│   ├── spain_tnbc
│   ├── codex_colon
│   ├── vectra_colon
│   ├── vectra_pancreas
│   ├── mibi_decidua
│   ├── mibi_breast
│   ├── experimental_results.csv
│   ├── gt_pred_ie_consolidated.csv
```
The `DCIS` and `spain_tnbc` folders and the `experimental_results.csv` file are part of this repository, the other folders (`codex_colon`, `vectra_colon`, `vectra_pancreas`, `mibi_decidua`) and  the `gt_pred_ie_consolidated.csv` file are available at [huggingface.co/datasets/JLrumberger/Pan-Multiplex](https://huggingface.co/datasets/JLrumberger/Pan-Multiplex). The folder `mibi_breast` will be made available through the huggingface data repository in the coming weeks.

### Dependencies
- Python 3.10
- Jupyter
- Pandas
- Numpy
- Scipy
- Matplotlib
- Seaborn
- Scikit-learn
- Scikit-image
- tqdm
- imageio
