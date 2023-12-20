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
Please see the libgen_notebook for detailed usage

#### Standard list preparation
The standard list mainly composes of 3 columns: name, inchikey, and mix label (the file names of the spectra).
Please refer to the sample standard list in the sample_data folder, and please do not change the column name.
With your sample list, please run following code to execute the preparation process. 

The prepared standard list will also show in standard_list directory, but with additional tail of '_cleaned' in file name

#### Feature finding
The feature finding is done by the custom code, ff_droup. It is a specialized feature finding process starting from MS/MS, rather than MS1.
The names of the raw spectra should be similar to the mix labels in the standard list, but does not to be identical since fuzzy match is enabled.
Also note that the files needs to be centroid, mzML files. If not, please use MSConvert.
Please run following code for the feature finding process.

#### Library curation
Major steps are match making (by precursor_mz), dereplication (based on ms1 intensity and SNR), recalibration (optional, but good to have if you have bad confidence in your mass spec), and spectral denoising.

#### Library exportation
Export library in msp or csv, whichever works for you the best!
Please also note that mgf and other popular formats are also supported, see detail in toolsets/file_io.py

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
