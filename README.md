# LibGen 2.0 stable version

This is the stable version of LibGen, with sample codes/data for showing the capability of the pipeline.

## Description

LibGen is a dedicated software for building high-quality mass spectral libraries from reference standard mixes.

## Getting Started

### Dependencies

The code was developed and Python 3.8. Run it in a dedicated conda environment is highly recommended. A yaml file will be uploaded soon.

### Installing

* The current version is only available by cloning the repository from GitHub.
* We plan to upload it onto PIP in the near future.


### Executing program
There are three main moudules in this pipeline: standard list preparation, feature finding and library curation.
* Creating a workspace
A workspace is just a folder for all of the files you will need for your library curation process. 
You can run the following code to create a workspace and put corresponding files in:
```
python setup_workspace.py [the_director_of_your_workspace]
```
The script will create a new folder with following subdirectories:
1. standard_list - the standard list you provided, as well as curated standard list by Standard list preparation section
2. spectra - all raw spectra files acquried from instrument
3. features - peak picking results resulted from feature finding section
4. curated library - the curated libraries
5. figures - saved figures (if any)

It will also print out the full directory to your workspace. Please use it in the following section.
* Standard list preparation
The standard list mainly composes of 3 columns: name, inchikey, and mix label (the file names of the spectra).
Please refer to the sample standard list in the sample_data folder, and please do not change the column name.
With your sample list, please run following code to execute the preparation process. 
```
python standard_list_prep.py [the_director_of_your_workspace]
```
The prepared standard list will also show in standard_list directory, but with additional tail of '_cleaned' in file name

* Feature finding
The feature finding is done by the custom code, ff_droup. It is a specialized feature finding process starting from MS/MS, rather than MS1.
The names of the raw spectra should be similar to the mix labels in the standard list, but does not to be identical since fuzzy match is enabled.
Please run following code for the feature finding process.
```
python feature_finding.py [the_director_of_your_workspace]
```
* Library curation
Now is eventually the time for curating the libarires! Please run the following code
```
python feature_finding.py -msp [the_director_of_your_workspace]
```
It will automatically export the curated libraries in both .msp format and .csv format. If you don't like .msp format, you can use the change -msp to:
1. -mgf: export to .mgf format
2. -ms: export to .ms file (SIRIUS propritery)
3. -mat: export to .mat file (MS-Finder propritery)

## Help

Please contact me directly for help information.

## Authors

Fanzhou Kong
Email: fzkong@ucdavis.edu

## Version History

* 1.0
    * Initial implementation of the project, deprecated.
* 2.0
    * Current release version of the LibGen, fully functioning.

## License

This project is licensed under Apache 2.0 license.

## Acknowledgments

Special thanks to Dr. Yuanyue Li, who developed the entropy similarity/spectral similarity algorithms, which serves as an essential part in this project.

Also appreciation goes to everyone in the Fiehn lab at UC Davis West Coast Metabolomics Center. This couldn't been doen without you.
