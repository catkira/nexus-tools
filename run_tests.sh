cd ..
rm -rf nexus-tools_build
mkdir nexus-tools_build
cd nexus-tools_build
export CXX=/usr/bin/g++-5
cmake -DLINUX_STATIC:BOOL=OFF -DCMAKE_BUILD_TYPE=Release ../nexus-tools
make test_adapter
make test_demultiplex
make test_general_processing
make test_trimming
./bin/test_adapter
./bin/test_demultiplex
./bin/test_general_processing
./bin/test_trimming

