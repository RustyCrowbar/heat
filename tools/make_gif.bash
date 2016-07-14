#!/usr/bin/env bash

# make_gif - Plot matrices graphs and make an animate gif out
# of it.
#
# Please run this script from the project's root directory!

typeset -lr img_dir="png_images"
typeset -la data_files=( $(ls | grep "rectangle_u") )
typeset -li i=${#data_files[@]}

mkdir -p "$img_dir"

echo "Plotting matrices..."
(( i = i - 1 ))
while (( i >= 0 ));
do
	local_file=\'${data_files[$i]}\'
	local_dir=\'$img_dir\'
	gnuplot -e "output_file_suffix=$i; output_dir=$local_dir; input_file=$local_file" \
		-c tools/plot_heatmap.plt
	(( i = i - 1 ))
done

output_files=$img_dir/*
convert -delay 30 -loop 0 $output_files $img_dir/heatmap.gif
