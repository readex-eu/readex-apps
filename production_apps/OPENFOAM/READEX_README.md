# 1] Repository content
 - OpenFOAM - OpenFOAM scripts and patch file
 - motorbikeXX - copy of a selected OpenFOAM test case


# 2] How to compile OpenFOAM
 * do not clone this repo to $HOME, because its size and number of files almost reaches the user cluster limits
 * download OpenFOAM-v1612+ with its ThirdParty libraries to OpenFOAM directory and apply the patch located in the directory (patch -p1 < $PATCH_FILE ): http://www.openfoam.com/download/release-history.php
 * compile the application using compile_1st.sh before you compile for other tools from the READEX workflow

