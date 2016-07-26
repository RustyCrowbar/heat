# plot.plt: Plot performance graphs using .dat files
#   given in parameter below.

set terminal png

outfile = "graph.png"
print outfile
set output outfile

set title 'Title'
set xlabel 'Matrix dimension'
set ylabel 'Time (s)'
# The key is where the curves' titles will be displayed
set key on
set key at 300,27
plot "seq.dat" using 1:2 title 'Sequential' with lines, \
     "parallel.dat" using 1:2 title 'Parallel' with lines, \
