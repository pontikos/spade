spade
=====

Modification of https://github.com/nolanlab/spade which works on a RData matrix instead of an FCS file.
The data matrices are stored in .RData files which are loaded.
Only one .RData file is required which is updated with a density (numeric), downsample (boolean) and cluster (factor) column, instead of storing these 3 separate files as was done in the original implementation.


