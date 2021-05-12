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
cd rpbasicdesign
python -m rpbasicdesign.cli --rpsbml_file tests/input/rp_1_12.sbml.xml
```

Output folders for dnabot-ready files and SBOL export can be set using `o_dnabot_dir` and `o_sbol_dir` options:
```sh
python -m rpbasicdesign.cli \
  --rpsbml_file tests/input/rp_1_12.sbml.xml \
  --o_dnabot_dir out/dnabot_input \
  --o_sbol_dir out/sbol_export
```

The number of constructs to design is tuned using `sample_size`:
```sh
python -m rpbasicdesign.cli \
  --rpsbml_file tests/input/rp_1_12.sbml.xml \
  --sample_size 5
```

Polycistronic constructs are built by default. To swtich to monocistronic, one can use the `monocistronic` option:
```sh
python -m rpbasicdesign.cli \
  --rpsbml_file tests/input/rp_1_12.sbml.xml \
  --monocistronic true
```

The complete list options is provided the embedded help, which can be printed using the `--help` or `-h` keywords:
```
python -m rpbasicdesign.cli -h

>>> usage: python -m dnaprep.cli [-h] --rpsbml_file RPSBML_FILE [--linker_parts_file LINKER_PARTS_FILE]
>>>                              [--linker_plate_file LINKER_PLATE_FILE] [--user_parts_file USER_PARTS_FILE]
>>>                              [--monocistronic MONOCISTRONIC] [--lms_id LMS_ID]
>>>                              [--lmp_id LMP_ID] [--backbone_id BACKBONE_ID] [--sample_size SAMPLE_SIZE]
>>>                              [--o_dnabot_dir O_DNABOT_DIR] [--o_sbol_dir O_SBOL_DIR]
>>>
>>> Convert rpSBML enzyme info in to BASIC construct. UniProt IDs corresponding enzyme variants are extracted from
>>> rpSBMl files. Promoters and RBSs are randomly chosen from a default list. CDSs, in other words gene
>>> variants, of enzymes are randomly chosen from amongst the UniProt IDs extracted. Constructs generated can be
>>> stored as (i) a CSV file ready to be used by DNA-Bot, (ii) as SBOL files.
>>>
>>> optional arguments:
>>>   -h, --help            show this help message and exit
>>>   --rpsbml_file RPSBML_FILE
>>>                         rpSBML file from which enzymes UniProt IDs will be collected
>>>   --linker_parts_file LINKER_PARTS_FILE
>>>                         File listing available linkers for constructs.
>>>   --linker_plate_file LINKER_PLATE_FILE
>>>                         File providing half linkers coordinates.
>>>   --user_parts_file USER_PARTS_FILE
>>>                         File listing user parts (eg backbone, promoters) available for constructs.
>>>   --monocistronic MONOCISTRONIC
>>>                         Build monocistronic constructs. Default to false, ie polycistronic constructs will
>>>                         be generated.
>>>   --lms_id LMS_ID       part ID to be used as the LMS methylated linker
>>>   --lmp_id LMP_ID       part ID to be used as the LMP methylated linker
>>>   --backbone_id BACKBONE_ID
>>>                         part ID to be used as the backbone
>>>   --sample_size SAMPLE_SIZE
>>>                         Number of construct to generate.
>>>   --o_dnabot_dir O_DNABOT_DIR
>>>                         Output folder to write construct and coord part files. It will be created if it does
>>>                         not exist yet. Existing files will be overwritten.
>>>   --o_sbol_dir O_SBOL_DIR
>>>                         Output folder to write SBOL depictions of constructs. It will be created if it does
>>>                         not exist yet. Existing files will be overwritten.
```

## Inputs

- `rpsbml_file` [MANDATORY]: SBML with retropath-like annotations. UnitProt IDs of enzyme are expected to 
  be listed here. More information of rpSBML file at [https://github.com/brsynth/rptools](https://github.com/brsynth/rptools).
  Some examples or rpSBML files are provided in `tests/input`.
- `linker_parts_file` [OPTIONAL]: CSV file listing the linker IDs available for the constructs. The content
  corresponding to the BioLegio commercial ([link](https://www.biolegio.com/products-services/basic/)) plate is
  used by default from the `rpbasicdesign/data/biolegio_parts.csv` file. A second predefined file corresponding to
  older version of the BioLegio plate is also described in `rpbasicdesign/data/legacy_parts.csv`.
- `linker_plate_file` [OPTIONAL]: CSV file providing the coordinates of half linker in the plate. By default, the
  `rpbasicdesign/data/biolegio_plate.csv` CSV file is used. An alternative file corresponding to an older version
  of the BioLegio plate is available in `rpbasicdesign/data/legacy_plate.csv`.
- `user_parts_file` [OPTIONAL]: CSV file listing user parts (eg backbone, promoters) available for constructs. The 
  content by default is ``
  
## For developers

### Installation from source code

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

### Predefined set of linkers

By default, the set of linkers used is the one presented available in from the commercial plate from BioLegio.
If one wants to use its own set of linkers, the user is advised to do it carefully and to look for more information.

### Linker naming conventions

Due to DNA-Bot implementation:
- RBS linkers should start with the `UTRn` suffix, where `n` could be any alphanumeric character.
- Any linkers should have its two half linkers ending with the `-P` and `-S` suffixes listed in the "plate" file, ie
  in the file that provides the well locations containing the DNA fragment. See the BASIC approach paper, and
  especially the supplementary files for more information.  

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
- Linkers have to be listed in `[biolegio|legacy]_parts.csv`.  
- Linker prefixes and suffixes well coordinates have to be listed in `[biolegio|legacy]_plate.csv`.
- Linkers described `user_parts.csv` are not considered.
- RBS linker IDs have to be in the form `AAA-BBB` with `AAA` being the linker suffix ID.

### Maximum number of CDSs per construct

The maximum number of genes in a construct limited to 3 with the default `biolegio_plate.csv` RBS library, because
there is only 3 different RBS suffix in the commercial BioLegio library.  

## TODO

- Better handle logs and add `verbose` option

## References

- Galaxy-SynBioCAD: https://doi.org/10.1101/2020.06.14.145730
- DNA-Bot: https://doi.org/10.1093/synbio/ysaa010
- BASIC assembly method: https://doi.org/10.1021/sb500356d
