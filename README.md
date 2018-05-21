# YSE Simulation Tools
---

### Conda environment

To make sure everyone is using the same version of Python,
there is an `environment.yml` that defines the content of a virtual
environment called `plasticcval`.

To create it first, run (only once)
```bash
conda env create -f environment.yml
```
Now everytime you need to use it
```bash
source activate YSE
```
You should be running Python 3.6 and everything will be fine.

If you need to leave the environment
```
source deactivate
```

### Download recent simulations
See [recent pkl file with simulations](https://www.dropbox.com/s/7mzf6xvmlwnbqyr/yse_ztf_YOUNG.pkl.gz?dl=0)

### Examples

Get histograms of the rise times of young simulated SNe
```python FindYoungSN.py --simname yse_ztf --empirical```

Get histograms of the rise times of just young SNe Ia
```python FindYoungSN.py --simname yse_ztf --empirical --sntype Ia```
