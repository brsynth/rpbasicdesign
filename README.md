# dnaprep

**A command-line tool to convert rpSBML files into SBOL and CSV files ready-to-be used with DNA-Bot".**
*dnaprep* extracts enzyme IDs from rpSBML files -- produced by the RP suite available in the [SynBioCAD Galaxy platform](https://galaxy-synbiocad.org) -- to generate genetic constructs compliant with the [BASIC](https://doi.org/10.1021/sb500356d) assembly approach.
CSV files produced then be used with [DNA-Bot](https://github.com/BASIC-DNA-ASSEMBLY/DNA-BOT) to generates python scripts allowing automated build of the genetic constructs using [OpenTrons](https://opentrons.com/) liquid handling robots.  
Formal description of parts involved in constructs is outputted in SBOL files.

## Install

From source code
```bash
git clone https://forgemia.inra.fr/thomas.duigou/dnaprep.git
cd dnaprep
conda env create -f environment.yaml -n <myenv>
pip install -e .
```

`<myenv>` has to be replaced by whatever meaningful name that will pleased the user.
If the no environment name is specified using the `-n` argument, the created environment will be `dnaprep`. 

## Usage

```bash
conda activate <myenv>
python -m dnaprep.cli --rpsbml_file tests/input/rp_1_12.sbml.xml --sample_size 12 --o_dnabot_file lala.csv --o_sbol_dir lala_sbol
```

Argument usage are described within the tool
```
python -m dnaprep.cli --help
```

## TODO

- Write tests

## References

- Galaxy-SynBioCAD: https://doi.org/10.1101/2020.06.14.145730
- DNA-Bot: https://doi.org/10.1093/synbio/ysaa010
- BASIC assembly method: https://doi.org/10.1021/sb500356d
 
