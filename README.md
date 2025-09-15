# RNA Fingerprinting Data

This package provides precomputed dictionaries for use with RNA fingerprinting. To install this package, you can run:

```
# install.packages("remotes"")
remotes::install_url("https://github.com/satijalab/rna-fingerprinting-data/raw/refs/heads/main/RNAFingerprintingData_0.1.0.tar.gz")
```

To view the available dictionaries, you can run the command `AvailableDictionaries()`. Then any desired dictionary can be installed with `InstallPrecomputedDictionary()`, and loaded in with `LoadPrecomputedDictionary`. For example:

```
library(RNAFingerprintingData)
AvailableDictionaries()
InstallPrecomputedDictionary(name = 'gwps')
gwps_dictionary <- LoadPrecomputedDictionary(name = 'gwps')
```