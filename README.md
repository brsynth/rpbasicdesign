# rpbasicdesign


[![Anaconda-Server Badge](https://anaconda.org/brsynth/rpbasicdesign/badges/latest_release_date.svg)](https://anaconda.org/brsynth/rxn_rebuild)
[![Anaconda-Server Badge](https://anaconda.org/brsynth/rpbasicdesign/badges/version.svg)](https://anaconda.org/brsynth/rxn_rebuild)

**A command-line tool to convert rpSBML files into SBOL and CSV files ready-to-be used with DNA-Bot.**

**rpbasicdesign** extracts enzyme IDs from rpSBML files -- produced by the RP suite available in the 
[SynBioCAD Galaxy platform](https://galaxy-synbiocad.org) -- to generate genetic constructs compliant with 
the [BASIC](https://doi.org/10.1021/sb500356d) assembly approach. CSV files produced are ready to be used with
[DNA-Bot](https://github.com/BASIC-DNA-ASSEMBLY/DNA-BOT) to generate instructions for automated build of the 
genetic constructs using [OpenTrons](https://opentrons.com/) liquid handling robots.

## Installation

```sh
conda install -c brsynth -c conda-forge rpbasicdesign
``` 


## Usage

Simple call:
```sh
conda activate <myenv>
python -m rpbasicdesign.cli --rpsbml_file tests/input/muconate_example.xml
```

Output folders for dnabot-ready files and SBOL export can be set using `o_dnabot_dir` and `o_sbol_dir` options:
```sh
python -m rpbasicdesign.cli \
  --rpsbml_file tests/input/muconate_example.xml \
  --o_dnabot_dir out/dnabot_input \
  --o_sbol_dir out/sbol_export
```

The number of constructs to design is tuned using `sample_size`:
```sh
python -m rpbasicdesign.cli \
  --rpsbml_file tests/input/muconate_example.xml \
  --sample_size 5
```

The complete list options is provided the embedded help, which can be printed using the `--help` or `-h` keywords:
```
python -m rpbasicdesign.cli -h

usage: python -m dnaprep.cli [-h]
                            --rpsbml_file RPSBML_FILE 
                            [--parts_files PARTS_FILES [PARTS_FILES ...]]
                            [--lms_id LMS_ID]
                            [--lmp_id LMP_ID]
                            [--backbone_id BACKBONE_ID]
                            [--sample_size SAMPLE_SIZE]
                            [--cds_permutation CDS_PERMUTATION]
                            [--o_dnabot_dir O_DNABOT_DIR]
                            [--o_sbol_dir O_SBOL_DIR]

Convert rpSBML enzyme info in to BASIC construct. UniProt IDs corresponding enzyme variants are extracted rpSBMl files. Promoters and RBSs are randomly chosen from a default list. CDSs, in other words gene variants,
of enzymes are randomly chosen from amongst the UniProt IDs extracted. Constructs generated can be stored as (i) a CSV file ready to be used by DNA-Bot, (ii) as SBOL files.

optional arguments:
  -h, --help            show this help message and exit
  --rpsbml_file RPSBML_FILE
                        rpSBML file from which enzymes UniProt IDs will be collected.
  --parts_files PARTS_FILES [PARTS_FILES ...]
                        List of files providing available linkers and user parts (backbone, promoters, ...) for constructs. Default: [data/biolegio_parts.csv, user_parts.csv]
  --lms_id LMS_ID       part ID to be used as the LMS methylated linker. Default: LMS
  --lmp_id LMP_ID       part ID to be used as the LMP methylated linker. Default: LMP
  --backbone_id BACKBONE_ID
                        part ID to be used as the backbone. Default: BASIC_SEVA_37_CmR-p15A.1
  --sample_size SAMPLE_SIZE
                        Number of construct to generate.Default: 3
  --cds_permutation CDS_PERMUTATION
                        Whether all combinations of CDS permutation should be built Default: true
  --o_dnabot_dir O_DNABOT_DIR
                        Output folder to write construct and plate files. It will be created if it does not exist yet. Existing files will be overwritten. Default: out/dnabot_in
  --o_sbol_dir O_SBOL_DIR
                        Output folder to write SBOL depictions of constructs. It will be created if it does not exist yet. Existing files will be overwritten. Default: out/sbol_export
```

## Inputs

This section documents input files required / optional, their purpose, and how information should be structured.

### rpSBML file [required]

SBML with retropath-like annotations. UnitProt IDs of enzyme are expected to be listed here. More information of rpSBML file at [https://github.com/brsynth/rptools](https://github.com/brsynth/rptools). Some examples or rpSBML files are provided in `tests/input`.

### Parts files [optional]

These are CSV files listing the linker IDs available for the constructs (BASIC linkers), as well as the user parts (backbone, promoters, ...). The format should be comma separated on 4 columns with header. Example below:
```
id,type,sequence,comment
L1,neutral linker,,
L2,neutral linker,,
L3,neutral linker,,
```

By default, the `rpbasicdesign/data/biolegio_parts.csv` file is used which corresponds to the BioLegio commercial plate ([link](https://www.biolegio.com/products-services/basic/)). A second predefined file corresponding to older version of the BioLegio plate is also described in `rpbasicdesign/data/legacy_parts.csv`.

For linkers, the `type` annotation should be one of `neutral linker`, `methylated linker`, `peptide fusion linker` or `RBS linker`. For user parts, `type` should be one of `backbone` or `constitutive promoter`. Other type will raise a warning and will be omited. By default, [biolegio_parts.csv](rpbasicdesign/data/biolegio_parts.csv) and [user_parts.csv](rpbasicdesign/data/user_parts.csv) are used.

Use the `parts_files` arguments to override.

**Important**:
- IDs should match the linker naming conventions (see below).
- IDs should match the IDs used in the plate file inputed to dnabot. As example -- but also ready to be used -- the [biolegio_plate.csv](rpbasicdesign/data/biolegio_plate.csv) is a valid input files for dnabot, with consistent IDs between `biolegio_parts.csv` and `biolegio_plate.csv`.

## Example

If one wishes to only use a subset of BASIC parts, the way to go is to provide a restricted list of parts with the `--parts_file` option.

The `examples` folder contains such an example:
```sh
python -m rpbasicdesign.cli \
    --rpsbml_file lycopene_example.xml \
    --sample_size 96 \
    --parts_file biolegio_parts_lycopene.csv \
    --o_dnabot_dir out/lycopene_dnabot \
    --o_sbol_dir out/lycopene_sbol
```

## For developers

### Installation

```sh
git clone https://github.com/brsynth/rpbasicdesign.git
cd rpbasicdesign
conda env create -f environment.yaml -n <myenv>
conda develop -n <myenv> .
```

### Tests

```sh
python -m pytest -v --cov=rpbasicdesign --cov-report html
```



## Constraints and limitations

The BASIC linker set is a major piece of the BASIC assembly method. For a detailed explanation of the BASIS approach,
see Storch et. al., ACS Synth. Biol., 2015 (doi: [10.1021/sb500356d](https://doi.org/10.1021/sb500356d)).

### Polycistronic constructs

Only polycistronic constructs are enabled at the moment.

### Predefined set of linkers

By default, the set of linkers used is the one presented available in from the commercial plate from BioLegio.
If one wants to use its own set of linkers, the user is advised to do it carefully and to look for more information.

### Linker naming conventions

Due to DNA-Bot implementation:
- RBS linkers should start with the `Un` suffix, where `n` could be any alphanumeric character.
- Any linkers should have its two half linkers ending with the `-P` and `-S` suffixes listed in the "plate" file, ie in the file that provides the well locations containing the DNA fragment. See the BASIC approach paper, and especially the supplementary files for more information.  

### Controlled vocabulary for parts file

Parts and linkers provided in the `*_parts.csv` files have to match on the following type:

- `neutral linker`
- `methylated linker`
- `RBS linker`
- `peptide fusion linker`
- `backbone`
- `constitutive promoter`

### Providing CDSs 

As of today, CDS are obtained only by parsing rpSBML files. 

### Custom linkers
For advanced users wishing to play with custom linkers:
- Linkers and parts can be provided using a custom file with the `--parts_files` argument.
- Linkers described `user_parts.csv` are not considered.
- RBS linker IDs have to be in the form `AAA-BBB` with `AAA` being the linker suffix ID.
- Linker prefixes and suffixes coordinates on the plate have to be listed in `[biolegio|legacy]_plate.csv`.

### Maximum number of CDSs per construct

The maximum number of genes in a construct limited to 3 with the default `biolegio_plate.csv` RBS library, because
there is only 3 different RBS suffix in the commercial BioLegio library.  

## TODO

- Better handle logs and add `verbose` option

## References

- Galaxy-SynBioCAD: https://doi.org/10.1101/2020.06.14.145730
- DNA-Bot: https://doi.org/10.1093/synbio/ysaa010
- BASIC assembly method: https://doi.org/10.1021/sb500356d
