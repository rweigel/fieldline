from setuptools import setup, find_packages

#install_requires = ["numpy","scipy","vtk","matplotlib", "urllib3"]
#            "git+https://github.com/GaryQ-physics/swmf_file_reader.git"]


install_requires = ["numpy","scipy","vtk","matplotlib", "urllib3",
                    "swmfio @ git+https://github.com/GaryQ-physics/swmfio.git"]

setup(
    name='fieldline',
    version='0.0.2',
    author='Gary Quaresima, Bob Weigel',
    author_email='rweigel@gmu.edu',
    packages=find_packages(),
    description='Field line tracing using multiple methods and libraries',
    install_requires=install_requires
   )
