from setuptools import setup, find_packages
import os
import re

_readme = 'README.md'
_extras = 'extras'

with open(_readme, 'r', encoding='utf-8') as f:
    _long_description = f.read()

with open(os.path.join(_extras, '.env'), 'r', encoding='utf-8') as f:
    for line in f:
        if line.startswith('PACKAGE='):
            _package = line.splitlines()[0].split('=')[1].lower()
        if line.startswith('URL='):
            _url = line.splitlines()[0].split('=')[1].lower()
        if line.startswith('AUTHORS='):
            _authors = line.splitlines()[0].split('=')[1].lower()
        if line.startswith('DESCR='):
            _descr = line.splitlines()[0].split('=')[1].lower()
        if line.startswith('CORR_AUTHOR='):
            _corr_author = line.splitlines()[0].split('=')[1].lower()

with open(os.path.join(_package, '_version.py'), 'r') as ifh:
    for line in ifh:
        m = re.search('__version__.*=.*"(.+)"', line)
        if m:
            _version = m.group(1)
            break


setup(
    name=_package,
    version=_version,
    author=_authors,
    author_email=_corr_author,
    description=_descr,
    long_description=_long_description,
    long_description_content_type='text/markdown',
    url=_url,
    packages=find_packages(),
    include_package_data=True,
    test_suite='pytest',
    license='MIT',
    classifiers=[
        'Topic :: Scientific/Engineering',
    ],
    keywords=[_package],
    # python_requires='>=3.6',
)
