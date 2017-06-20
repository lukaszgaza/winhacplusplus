# Winhac++

### Preparing Linux CERN 6 (SLC6) environment
1. Install xerces-c
```
sudo yum install xerces-c
sudo yum install xerces-c-devel
```

2. Install LHAPDF
```
sudo yum install lhapdf
sudo yum install lhapdf-devel
sudo yum install lhapdf-pdfsets-minimal
```

3. Install Boost
```
sudo yum install boost
sudo yum install boost-devel
```

4. Install HepMC
```
sudo yum install HepMC
sudo yum install HepMC-devel
```

5. Install Pythia8
```
sudo yum install pythia8
sudo yum install pythia8-devel
```

6. Install Root
```
sudo yum install root
sudo yum install root-roofit
sudo yum install root-physics
```

7. Install development tools
```
sudo yum install cmake
sudo yum install gcc-g++
sudo yum install gcc-gfortan
```

### Building under Linux CERN 6 (SLC6)
```
cmake .
make
```
