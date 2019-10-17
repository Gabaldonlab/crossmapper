# To Build Crossmapper Package
```bash
## run in working directory path/to/build/
## for python 3.6
conda build crossmapper -c bioconda --python=3.6
## for python 3.7
conda build crossmapper -c bioconda --python=3.7

## then convert to osx
conda convert -p osx-64 path/to/package*

## to upload to conda channel do for all packages
anaconda upload path/to/package*
```
