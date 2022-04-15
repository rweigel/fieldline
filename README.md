# Overview

Python package for tracing field lines using multiple methods and libraries

Tested for python 2.7.18 and 3.7.9

# Install

This package has an optional dependency on "swmf_file_reader", which will not be automatically installed since it is not on PyPi.

As such, you need to follow the installation instructions at [https://github.com/GaryQ-physics/swmf_file_reader](https://github.com/GaryQ-physics/swmf_file_reader) before continuing.

`swmf_file_reader` is used to convert SWMF's output files to VTK files for fieldline tracing in ParaView. 

# Instructions

User
```
pip install 'git+https://github.com/rweigel/fieldline' --upgrade
python -c 'from fieldline import demos; demos.demo2()
```

Developer
```
git clone https://github.com/rweigel/fieldline
cd fieldline
pip install --editable .
```
