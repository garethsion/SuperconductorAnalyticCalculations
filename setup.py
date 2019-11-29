from setuptools import setup, find_packages
#from setuptools.command.develop import develop
#from setuptools.command.install import install

setup(
        name="supercalcs",
        version="0.1.0",
        packages=find_packages(exclude=['*test','*ipynb']),
        author="Gareth Sion Jones",
        author_email="garethsion@googlemail.com",
        python_requires='>3.4.0',
        install_requires=['numpy','scipy', 'pandas', 'scikit-rf','matplotlib']
    ) 
