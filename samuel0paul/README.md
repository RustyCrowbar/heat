2D-Heat-Transfer
================
##Theory - Physics
When we divide a 2-D plate which is in steady state condition(ie. there is no more temperature variation in a point with respect to time), We find:

```Matlab
temp = f(posX, posY)
% _approx: ~=
temp(m, n) ~= (temp(m, n-1) + temp(m, n+1) + temp(m-1, n) + temp(m+1, n)) / 4
```

The above postulate stands true if the heat transfer is only through Conduction(i.e. there is no heat transfer through convection or radiation).

Implementation
--------------

What has been done to try to enhance performances (in chronological order):
* *Use TBB _parallel for_* to parallelize the computation
  for inner nodes only. Outer nodes are left in the serial
  world because taking them into account would require a lot of
  testing (_if_) and break the processor pipeline.
  This works pretty well as a starting point.
* *Make it possible to compute multiple iterations* within a single
  _parallel for_. This was abandonned because it required taking
  account of outer nodes, and filling the _parallel for_ with
  branch conditions was not considered ideal.
  We tried to use 2 _std::threads_ and make them communicate:
  1 would compute for outer nodes and copying, the other for
  inner nodes. The problem is, we would need to give TBB our own
  chunking algorithm so the other thread would know how synchronize
  properly. This is possible, but is quite ugly.
* *Get tricky and actually use a larger matrix* to put everything
  in the _parallel for_ without annoying branch conditions.
  This was tried by using 3-uples instead of pairs, the 3rd element
  being the number of neighboring nodes. All extra outer nodes would
  be zeroed-valued, thus not impacting the calculation.
  The resulting code was put in the _parallel fail_ branch, and
  did not prove efficient.
* *Put every node calculation in the _parallel for_*, and try to
  implement easily multiple-iterations-in-one-go. The problem was:
  how to synchonize the different subranges that TBB generated?
  In CUDA for example, we could have used _syncthreads()_.
  Here, TBB's _parallel for_ does not provide such synchronization
  mechanism. Thus we stopped working this way. The incomplete code
  is in the _big parallel for_ branch.
* *Rewrite the parallel version using TBB graphs.*
  More specifically, using a TBB dependance graph. This allows parallelism
  between calculations for inner and outer nodes. Even more important, it
  is easy to make to graph a loop and stop cleanly at any number of iteration.
  This allows charging up the CPU's threads and provides drastic
  performance increasments.

Various ideas on how to better use Intel TBB:

### Using a different partitioner (obsolete)

By default, the _automatic partitioner_ is used. We could set it to the
_simple partitioner_ and manually tweak the limit. We could also try
the _affinity partitioner_ but the same data are not being use more than
a few times.

### Changing the grain size (obsolete)

And plot graphs to see what changed (with different problem sizes).


LICENSE
-------

```license
The MIT License (MIT)

Copyright (c) 2014 Samuel Vishesh Paul

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
