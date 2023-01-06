# Release history

## 1.2.1 (2023-01-06)

- chore(test_Construct.py): compare based on Component URIs
- chore(test.yml): transition to mamba for GA
- chore: fix github action workflows

## 1.2.0 (2022-12-05)

- chore(Designer): add type hint
- fix(Designer): prevent StopIteration exception
- fix(Designer): follow csv writer end of line recommandation

## 1.1.1 (2022-08-23)

- fix(test_Construct.py): update hash probably due to a new version of lxml library
- test(test_get_sbol_new.txt): add an output test file for debugging (function test_get_sbol)

## 1.1.0 (2022-05-05)

- feat(cli): expose max number of gene to CLI args
- docs(README): update
- chore: remove deprecated code

## 1.0.1 (2022-03-25)

- fix(test_Construct.py): update hash for sbol output

## 1.0.0 (2022-03-23)

- feat(cli.py): SBOL output dir is now optional
- fix(cli.py): update default number of constructs
- docs(README): update
- build(meta.yaml): remove unneeded library

## 0.3.4 (2021-10-25)

- docs: add README about the lycopene file
- feat: get enzyme IDs from selenzy annotation
- chore: update example with new selenzy output file
- build: make use of rptools library

## 0.3.3 (2021-10-22)

- fix: handle correctly multiple UniProt IDs per reaction
- fix: ID based workaround to generate duplicate parts in SBOL

## 0.3.2 (2021-10-15)

- fix: generate plate coordinates by columns
- fix: output up to 88 constructs, not 96
- docs: provide example for the lycopene pathway

## 0.3.1 (2021-07-22)

- style(cli): show default values in help
- feat(cli): add cds_permutation arg
- feat(Designer): export the biolegio plate file

## 0.3.0 (2021-07-20)

- docs(README): update info on RBS naming pattern
- docs: provides an example of the --parts_file option
- fix(Designer): fix the way combinatorics is built
- !fix: constructs are always polycistronic

## 0.2.2

- fix(Designer): raise exception if no linkers / parts provided

## 0.2.1

- fix(Designer): _parts_files not defined

## 0.2.0

- feat: single argument for all part files (linker, user)
- feat: remove the need for the linker plate file
- fix(Designer): uniprot ID extraction
- docs(README): update input section
- docs(README): update command line example

## 0.1.1

- test(test_designer.py): fix unsuccessful test (output file renamed)
- build(environment.yaml): fix sbol2 package not found (due to package new name)
- fix(designer): fix spelling in output files
- build(recipe): update dependencies
- build: no longer need brsynth conda channel, pysbol2 is in conda-forge

## 0.1.0

- fix(designer): fix random generation not "reproductible" with a fixed seed
- refactor: separate classes
- fix: fix loading of default data
- docs(readme): add developer section
- fix: update dependancies to sbol2 as pysbol is deprecated
- refactor(cli): build polycistronic constructs by default
- feat: refine part and linker attributes, predict linker suffix and prefix

## 0.0.2

- Packaged version 

## 0.0.1

- Hello, world
