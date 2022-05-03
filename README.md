# kdeLF

**kdeLF** is an MIT licensed Python implementation of Yuan et al.’s [Kernel Density Estimation (KDE) method for estimating luminosity functions](https://arxiv.org/abs/2203.06700). It is a wrapper to compiled fortran code that does the heavy lifting, and is therefore relatively fast. We are open to all questions, feedback, commentary, and suggestions as long as they are constructive. Discussions should always come in the form of git issues. 

### Documentation

Read the docs at [kdelf.readthedocs.io](https://kdelf.readthedocs.io/en/latest/).

## Installation and test

### Using pip

The recommended way to install the stable version of *kdeLF* is using [pip](http://www.pip-installer.org/):

```
pip install -U kdeLF
```

We have uploaded several `.whl` files to the [PyPI](https://pypi.org/project/kdeLF) web to support as many platforms as possible. If your platforms happens to be an exception, then the pip installation may fail. In this situation, you need to install the [Intel fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/news/free-intel-software-developer-tools.html) first, and then try the pip installation again.  In addition, the version of numpy should be 1.21.0 If you have problems installing, please open an issue at [GitHub](https://github.com/yuanzunli/kdeLF/issues).

### From source

You can also install *kdeLF* after a download from [GitHub](https://github.com/yuanzunli/kdeLF/). Note that this requires a [Intel fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/news/free-intel-software-developer-tools.html) to be installed in advance.

```
git clone https://github.com/yuanzunli/kdeLF.git
cd kdeLF
pip install .
```

### Test the installation

To make sure that the installation went alright, you can execute the built-in test program in **kdeLF** by running the following command:

```
python3 -m kdeLF.test_kdeLF
```

## Citation

Please cite the following papers if you found this code useful in your research:
1. Yuan, Z., Zhang, X., Wang, J., Cheng, X., & Wang, W. 2022, ApJS, 248, 1 ([arXiv](https://arxiv.org/abs/2203.06700), [ADS](https://ui.adsabs.harvard.edu/abs/2020ApJS..248....1Y), [BibTeX](https://ui.adsabs.harvard.edu/abs/2020ApJS..248....1Y/exportcitation)).
2. Yuan, Z., Jarvis, M. J., & Wang, J. 2020, ApJS, 248, 1 ([arXiv](https://arxiv.org/abs/2003.13373), [ADS](https://ui.adsabs.harvard.edu/abs/2020ApJS..248....1Y), [BibTeX](https://ui.adsabs.harvard.edu/abs/2020ApJS..248....1Y/exportcitation)).


