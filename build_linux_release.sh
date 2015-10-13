cd ..
mkdir nexus-tools_build
cd nexus-tools_build
export CXX=/usr/bin/g++-5
cmake -DLINUX_STATIC:BOOL=ON -DCMAKE_BUILD_TYPE=Release ../nexus-tools
make flexbar++
make nexus-pre
make MappingAnalyzer
cp -f bin/flexbar++ ../nexus-tools/bin/
cp -f bin/nexus-pre ../nexus-tools/bin/
cp -f bin/MappingAnalyzer ../nexus-tools/bin/


