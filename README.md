# fieldline
Python package for tracing field lines using multiple methods and libraries

Tested for python 2.7.18 and 3.7.9

## install

### Dependency Note
this package has an optional dependency on "swmf_file_reader",
which will not be automatically installed since it is not on PyPi.

As such, you need to follow the installation instructions at
![https://github.com/GaryQ-physics/swmf_file_reader](https://github.com/GaryQ-physics/swmf_file_reader)
before continuing

"swmf_file_reader" is used to convert swmf's output files to vtk for vtk fieldline tracing.
If this functionality isn't desired, then it can be left out.

### Instructions
for user:
```
pip install 'git+https://github.com/rweigel/fieldline' --upgrade
python -c 'from fieldline import demos; demos.demo2()
```

for developer:
```
git clone https://github.com/rweigel/fieldline
cd fieldline
pip install --editable .
```
