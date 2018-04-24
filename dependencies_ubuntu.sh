sudo apt-get install -y libgsl-dev libfftw3-dev libpng++-dev libcairo2-dev freeglut3-dev libtiff5-dev qt4-default libqt4-dev libcgal-dev libgsl-dbg liblas-c-dev liblas-dev libblas-dev
git clone https://github.com/pngwriter/pngwriter.git
mkdir -p pngwriter-build
cd pngwriter-build
cmake ../pngwriter
make -j
sudo make install

 

