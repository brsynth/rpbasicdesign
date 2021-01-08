# rpbasicdesign

**A command-line tool to convert rpSBML files into SBOL and CSV files ready-to-be used with DNA-Bot.**
*rpbasicdesign* extracts enzyme IDs from rpSBML files -- produced by the RP suite available in the [SynBioCAD Galaxy platform](https://galaxy-synbiocad.org) -- to generate genetic constructs compliant with the [BASIC](https://doi.org/10.1021/sb500356d) assembly approach.
CSV files produced then be used with [DNA-Bot](https://github.com/BASIC-DNA-ASSEMBLY/DNA-BOT) to generates python scripts allowing automated build of the genetic constructs using [OpenTrons](https://opentrons.com/) liquid handling robots.  
Formal description of parts involved in constructs is outputted in SBOL files.

## Install

From source code
```bash
git clone https://github.com/brsynth/rpbasicdesign.git
cd rpbasicdesign
conda env create -f environment.yaml -n <myenv>
pip install -e .
```

`<myenv>` has to be replaced by whatever meaningful name that will pleased the user.
If the no environment name is specified using the `-n` argument, the created environment will be `rpbasicdesign`. 

## Usage

```bash
conda activate <myenv>
python -m rpbasicdesign.cli --rpsbml_file tests/input/rp_1_12.sbml.xml --sample_size 12 --o_dnabot_file lala.csv --o_sbol_dir lala_sbol
```

Argument usage are described within the tool
```
python -m rpbasicdesign.cli --help
```

## BASIC linkers

The BASIC linker set is a major piece of the BASIC assembly method. For a detailed explanation of the BASIS approach, see Storch et. al., ACS Synth. Biol., 2015 (doi: [10.1021/sb500356d](https://doi.org/10.1021/sb500356d)).

### Predefined set of linkers

By default, the set of linkers used is the one presented in the original BASIC approach.
If one wants to use its own set of linkers, the user is advised to do it carefully and to look for more information.

### Linker naming naming conventions

Due to DNA-Bot implementation:
- RBS linkers should start with the `UTRn` suffix, where could be any alpha numeric character.
- Any linkers should have its two half linkers ending with the `-P` and `-S` suffixes listed in the "part file", ie in the file that provides the well locations containing the DNA fragment. See the BASIC approach paper, and especially the supplementary files for more information.  

## Default parts

Default backbone, linkers and promoters are used. Definition of the parts are in `rpbasicdesign/data/default_parts.tsv`.

## TODO

- Write tests
- Better handle logs and add `verbose` option
- Conda packaging
- CI

## Constraints and limitations

### Controlled vocabulary

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

## References

- Galaxy-SynBioCAD: https://doi.org/10.1101/2020.06.14.145730
- DNA-Bot: https://doi.org/10.1093/synbio/ysaa010
- BASIC assembly method: https://doi.org/10.1021/sb500356d
 
