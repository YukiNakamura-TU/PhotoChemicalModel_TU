#!/bin/csh -f
#!/bin/csh -f
#PJM -L "elapse=144:00:00"
#PJM -L "vnode=1"
#PJM -L "vnode-core=8"
#PJM -j
#PJM -X

cd build
cmake ..
cmake --build .
cd ..
mpirun -np 8 ./PhotoChemistry