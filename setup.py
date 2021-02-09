from setuptools import setup, find_packages

install_requires = ["vtk"]

setup(
    name='fieldline',
    version='0.0.1',
    author='Gary Quaresima, Bob Weigel',
    author_email='rweigel@gmu.edu',
    packages=find_packages(),
    description='Field line tracing using multiple methods and libraries',
    install_requires=install_requires
