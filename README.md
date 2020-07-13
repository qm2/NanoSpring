
#### Compilation:
Put nvcc in PATH.

```
sudo apt-get libgoogle-perftools-dev
mkdir bin/
make -j
```

#### Use
```
python createData.py
./testMinHash test.reads k n # set k and n
# when this stops input overlapBaseThreshold numKMersThreshold
```

TODO: 
- make usable both with and without GPUs and CUDA
- make compilation install dependencies automatically (e.g., gperftools)
