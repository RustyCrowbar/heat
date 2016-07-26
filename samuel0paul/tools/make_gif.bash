#!/usr/bin/env bash

# make_gif - Plot matrices graphs and make an animate gif out
# of it. The script expects text files of the following format
# to be present at the current directory:
#   rectangle_u[0-9][0-9][0-9][0-9].txt
#
# Please run this script from the project's root directory!

set -e

typeset -lr img_dir="png_images"
typeset -la data_files=( $(ls | grep "temperature_*") )
typeset -li i=${#data_files[@]}

mkdir -p "$img_dir"
echo -n "Plotting matrices... "

cores=$(grep -c ^processor /proc/cpuinfo)
(( cores *= 2 ))

jobs=0

for (( i = i - 1; i >= 0; i = i - 1 )); do
	local_file=\'${data_files[$i]}\'
	local_dir=\'$img_dir\'


        while (( jobs > cores )); do
            jobs=$(jobs -r | wc -l)
        done

	gnuplot -e "output_file_suffix=$i; output_dir=$local_dir; input_file=$local_file" \
		-c tools/plot_heatmap.plt &

        (( ++jobs ))
done &>/dev/null

echo "OK"

while (( jobs > 0 )); do
    jobs=$(jobs -r | wc -l)
done

echo -n "Building animated GIF... "
output_files=$img_dir/*
convert -delay 15 -loop 0 $output_files heatmap.gif

echo "OK"

rm -rf $img_dir
