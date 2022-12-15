# fieldline

Python package for tracing field lines using multiple methods and libraries

Tested on Python 2.7.18 and 3.7.9

## Install

### Instructions

User:

```
pip install 'git+https://github.com/rweigel/fieldline' --upgrade
python -c 'from fieldline import demos; demos.demo2()
```

Developer:

```
git clone https://github.com/rweigel/fieldline
cd fieldline
pip install --editable .
```

### Dependency Note

This package has an optional dependency on "swmf_file_reader",
which will not be automatically installed since it is not on PyPi.

As such, follow the installation instructions at
![https://github.com/GaryQ-physics/swmf_file_reader](https://github.com/GaryQ-physics/swmf_file_reader)
before continuing.

"swmf_file_reader" is used to convert SWMF's output files to VTK for fieldline
tracing using the VTK library. If this functionality isn't desired, then it
can be left out.

## See also

# Overview

Python package for tracing field lines using multiple methods and libraries

