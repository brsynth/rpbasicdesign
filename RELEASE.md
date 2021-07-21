# Release history

## 0.3.1
- style(cli): show default values in help
- feat(cli): add cds_permutation arg
- feat(Designer): export the biolegio plate file

## 0.3.0
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


