#!usr/bin/bash 

# Coded by Jingjing liu (jingjing.liu@stjude.org)

# python: 3.6.1
# SJARACNe: https://github.com/jyyulab/SJARACNe

# This script runs on files in the current working directory

sjaracne local -e ./SJAR/input.exp -g ./SJAR/tf.txt -n 2 -o ./SJAR/outputs/cwl/cwltool/SJARACNE_out.final
sjaracne local -e ./SJAR/input.exp -g ./SJAR/sig.txt -n 2 -o ./SJAR/outputs/cwl/cwltool/SJARACNE_out.final
