# plot_heatmap.plt: Draw a matrix representing the heatmap described
#   in the input file and save it in a png file.
#
# output_file_suffix: Number (3 digits max) to append to the output filename
# output_dir        : Output directory for png files
# input_file        : Input .txt file to use for matrix plotting

set terminal png
outfile = sprintf('%s/heat_%d.png', output_dir, output_file_suffix)
print outfile
set output outfile
plot input_file matrix with image
