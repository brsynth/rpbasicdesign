# dnaprep

**A command-line tool to convert rpSBML files into SBOL and CSV files ready-to-be used with DNA-Bot"**
*dnaprep* extracts enzyme IDs from rpSBML files -- produced by the RP suite available in the [SynBioCAD Galaxy platform](https://galaxy-synbiocad.org) -- to generate genetic constructs compliant with the [BASIC](https://doi.org/10.1021/sb500356d) assembly approach.
CSV files produced then be used with [DNA-Bot](https://github.com/BASIC-DNA-ASSEMBLY/DNA-BOT) to generates python scripts allowing automated build of the genetic constructs using [OpenTrons](https://opentrons.com/) liquid handling robots.  
Formal description of parts involved in constructs is outputted in SBOL files.

## Install

```bash
conda create -n <myenv> python=3
conda activate <myenv>
conda install -c conda-forge rdkit=2020 jupyterlab python-libsbml
pip install brs_libs=0.3.1
```
`<myenv>` has to be replaced by whatever meaningful name that will pleased the user.

## Usage

```bash
conda activate <myenv>
  ..
```

## References

- Galaxy-SynBioCAD: https://doi.org/10.1101/2020.06.14.145730
- DNA-Bot: https://doi.org/10.1093/synbio/ysaa010
- BASIC assembly method: https://doi.org/10.1021/sb500356d
 
