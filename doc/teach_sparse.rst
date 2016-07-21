.. HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
.. HQ X
.. HQ X   quippy: Python interface to QUIP atomistic simulation library
.. HQ X
.. HQ X   Copyright James Kermode 2010
.. HQ X
.. HQ X   These portions of the source code are released under the GNU General
.. HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
.. HQ X
.. HQ X   If you would like to license the source code under different terms,
.. HQ X   please contact James Kermode, james.kermode@gmail.com
.. HQ X
.. HQ X   When using this software, please cite the following reference:
.. HQ X
.. HQ X   http://www.jrkermode.co.uk/quippy
.. HQ X
.. HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Command line options of the teach_sparse main program
=====================================================

Main options
------------
.. autofunction:: quippy.teach_sparse_parse_command_line

GAP options
-----------
.. autofunction:: quippy.teach_sparse_parse_gap_str

`sparse_method` options are:
 - RANDOM: default, chooses n_sparse random datapoints
 - PIVOT: based on the full covariance matrix finds the n_sparse "pivoting" points
 - CLUSTER: based on the full covariance matrix performs a k-medoid clustering into n_sparse clusters, returning the medoids
 - UNIFORM: makes a histogram of the data based on n_sparse and returns a data point from each bin
 - KMEANS: k-means clustering based on the data points
 - COVARIANCE: greedy data point selection based on the sparse covariance matrix, to minimise the GP variance of all datapoints
 - UNIQ: selects unique datapoints from the dataset
 - FUZZY: fuzzy k-means clustering
 - FILE: reads sparse points from a file
 - INDEX_FILE: reads indices of sparse points from a file
 - CUR_COVARIANCE: CUR, based on the full covariance matrix
 - CUR_POINTS: CUR, based on the datapoints
