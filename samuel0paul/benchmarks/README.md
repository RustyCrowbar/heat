The *.dat* files are data files (like .csv files basically)
for different types of parallelism level (sequential,
parallel with 1, 2, 4, or 10 granularity).

They were fed to _gnuplot_ to obtain the *.png* images
listed in each subdirectory.

Subdirectories:
* "change matrix dim" compares the different implementations
  on growing square matrix dimensions.
* "small matrix dim" does the same, but with very small (<=5O)
  dimensions.
* "old implem" compares the sequential and parallel forms,
  when the parallel one was just using a parallel for and
  a parallel reduce.
* "small old implem" does the same, but with very small (<=5O)
  dimensions.
 
