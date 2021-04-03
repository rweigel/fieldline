# fieldline
Python package for tracing field lines using multiple methods and libraries

note, currently writing vtk only works with python 2

## install

Note: for using vtk tracing, this package has dependency on "swmf_file_reader",
which will not be automatically installed since it is not on PyPi.

As such, you need to follow the installation instructions at
![https://github.com/GaryQ-physics/swmf_file_reader](https://github.com/GaryQ-physics/swmf_file_reader)
befor continuing

for user:
```
pip install 'git+https://github.com/rweigel/fieldline' --upgrade
cd fieldline
python trace_demo.py
```

for developer:
```
git clone https://github.com/rweigel/fieldline
cd fieldline
pip install --editable .
python trace_demo.py
```
