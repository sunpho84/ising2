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
OMP_NUM_THREADS=2 ./ising
```

