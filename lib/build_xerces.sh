tar -zxvf xerces-c-3.0.1.tar.gz
cd xerces-c-3.0.1
./configure --prefix=$PWD/build
make
make install