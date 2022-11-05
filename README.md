#  synthetic_migraine_data

This repo is made as proof of concept.


### class Syndata in get_syndata.py:

Synthetic Migrane-Days are simulated by adding colored noise to a 
biorythm. The biorythm itself is assumed to be sub-threshold, but additional
colored noise is triggering migraine resulting in migraine attack.

Synthetic data are binary data, where 1 stands for "headache day" and 0 stands for "non-headache day".


### plot_spectram.py:

uses Syndata class to create syntheic data.

Spectral analysis is performed in order to retrieve the frequency and frequency modulation of the
underlying biorythm
