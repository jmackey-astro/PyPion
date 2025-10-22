

# For development

## local compilation and installation

```
python3 -m venv ~/.local/venv
source ~/.local/venv/bin/activate
python3 -m pip install astropy
python3 -m pip install numpy pandas matplotlib 

# install pypion into the virtual env:
git clone git@github.com:greensh16/PyPion.git 
cd PyPion/
git checkout 3-add-yt-db-creator
#python3 setup.py build
#python3 setup.py install
python -m pip install ./
```

This worked for JM on macOS with homebrew python 3.13

