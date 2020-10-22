from setuptools import setup, find_packages

setup(
    name='dnaprep',
    version='0.0.1',
    description='Convert rpSBML enzyme info in to BASIC construct.',
    license='MIT',
    author='Thomas Duigou',
    author_email='thomas.duigou@inrae.fr',
    url='https://github.com/brsynth/dnaprep',
    packages=find_packages(),
    keywords=['dnaprep'],
    classifiers=[
        'Topic :: Scientific/Engineering',
    ]
)
