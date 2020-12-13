#!/bin/bash
#

# beam is $1

./beamformer -d ~/data/fl_corr01.out -f ~/data/beamformer_weights_corr01.dat -o ~/data/corr01.fil -z 1498.75 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr02.out -f ~/data/beamformer_weights_corr02.dat -o ~/data/corr02.fil -z 1487.03125 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr03.out -f ~/data/beamformer_weights_corr03.dat -o ~/data/corr03.fil -z 1475.3125 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr21.out -f ~/data/beamformer_weights_corr21.dat -o ~/data/corr21.fil -z 1463.59375 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr05.out -f ~/data/beamformer_weights_corr05.dat -o ~/data/corr05.fil -z 1451.875 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr06.out -f ~/data/beamformer_weights_corr06.dat -o ~/data/corr06.fil -z 1440.15625 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr07.out -f ~/data/beamformer_weights_corr07.dat -o ~/data/corr07.fil -z 1428.4375 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr08.out -f ~/data/beamformer_weights_corr08.dat -o ~/data/corr08.fil -z 1416.71875 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr09.out -f ~/data/beamformer_weights_corr09.dat -o ~/data/corr09.fil -z 1405.0 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr10.out -f ~/data/beamformer_weights_corr10.dat -o ~/data/corr10.fil -z 1393.28125 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr11.out -f ~/data/beamformer_weights_corr11.dat -o ~/data/corr11.fil -z 1381.5625 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr12.out -f ~/data/beamformer_weights_corr12.dat -o ~/data/corr12.fil -z 1369.84375 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr13.out -f ~/data/beamformer_weights_corr13.dat -o ~/data/corr13.fil -z 1358.125 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr14.out -f ~/data/beamformer_weights_corr14.dat -o ~/data/corr14.fil -z 1346.40625 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr15.out -f ~/data/beamformer_weights_corr15.dat -o ~/data/corr15.fil -z 1334.6875 -s 1.4 -n $1
./beamformer -d ~/data/fl_corr16.out -f ~/data/beamformer_weights_corr16.dat -o ~/data/corr16.fil -z 1322.96875 -s 1.4 -n $1

python splicer.py 55000.0 ${1} ~/data/spliced.fil ~/data/corr01.fil ~/data/corr02.fil ~/data/corr03.fil ~/data/corr21.fil ~/data/corr05.fil ~/data/corr06.fil ~/data/corr07.fil ~/data/corr08.fil ~/data/corr09.fil ~/data/corr10.fil ~/data/corr11.fil ~/data/corr12.fil ~/data/corr13.fil ~/data/corr14.fil ~/data/corr15.fil ~/data/corr16.fil 



