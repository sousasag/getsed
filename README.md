# getsed

## Get the code

git clone https://github.com/sousasag/getsed

## Update the path in the code (line 10) in getsed.py


## Get SED Kurucz Models

### sedfitter

```
wget ftp://ftp.astro.wisc.edu/outgoing/tom/model_packages/models_kurucz_05sep11.tgz
tar xvf models_kurucz_05sep11.tgz
```

### stsci ck04 models
```
wget -r ftp://ftp.stsci.edu/cdbs/grid/ck04models
mv ftp.stsci.edu/cdbs/grid/ck04models ../
```

## Using the code

Add the folder to your Python path.
Some examples:
```
export PYTHONPATH="${PYTHONPATH}:/my/path/to/getsed/"
```
or
```
PYTHONPATH=$PYTHONPATH:/my/path/to/getsed/
```
or inside a python script:
```
import sys
sys.path.append('/my/path/to/getsed/')
import getsed

wave, flux = get_sed_interpolated_cube(5777, 0.0, 4.4)
```

