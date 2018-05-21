from setuptools import setup, find_packages

with open("README.md") as flines:
    readme = flines.read()

with open("LICENSE") as flines:
    license = flines.read()

with open("requirements.txt") as flines:
    requirements = [line.strip() for line in flines] 

setup(
    name='proteinnetworks',
    version='0.1dev',
    description='Network analysis tools for protein structure',
    long_description=readme,
    url='http://github.com/wllgrnt/protein-networks',
    author='William Grant',
    author_email='contact@wpg.io',
    license=license,
    setup_requires=['pytest_runner'],
    tests_require=['pytest'],
    install_requires=requirements,
    packages=find_packages("src", exclude=('tests', 'examples','htmlcov')),
    package_dir={'': 'src'}
)
