Compiling
---

Bare compilation:

```c++
g++ -o ising ising.cpp 
```

Enable optimization:
```
g++ -o ising ising.cpp -O3
```

Enable plot (needs gnuplot):
```
g++ -o ising ising.cpp -DPLOT
```

Enable openmp:
```
g++ -o ising ising.cpp -DPLOT -fopenmp
```

Parallel execution:
```
OMP_NUM_THREADS=2 ./ising L nConfs
```

Assignement Lect 5
-------------------

1) description of the system
2) description of the algorithm
3) plot of the simulation in the cold/hot phase
4) scaling of the time with the volume in the range L = 8 - 128
5) parallelization of the energy calculation: verify the reproducibility (same result is obtained every time the code is run)
6) scaling of the time with the number of threads in two regimes: small and large volume, with different placement on cpu
7) point out the different behavior in the two regimes


OMP_verbosity and placement
-------------
```
export OMP_DISPLAY_AFFINITY=true
export OMP_DISPLAY_ENV=true
OMP_NUM_THREADS=4 OMP_PROC_BIND=true ./ising L nConfs
OMP_NUM_THREADS=4 OMP_PROC_BIND=true OMP_PLACES={0,48} ./ising L nConfs
```
in the list of places, the id of the core to be used is labelled with P#


CPU structure
---

![cpu](cpu.svg)


Running different parameters
---

```
$ cat pars.txt #L nconfs nthreads
3 100 1
3 100 2
3 100 3
3 100 4
...

$ while read L nconfs nthreads;do OMP_NUM_THREADS=$nthreads ./ising $L $nconfs;done < pars.txt
```

