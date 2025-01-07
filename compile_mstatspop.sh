#To compile:
brew install gsl
brew install htslib
brew install zlib
sh build.sh
cp ./build/mstatspop ./bin
cp ./build/tfa_merge ./bin
cp ./build/tfa_index ./bin
cp ./build/ms ./bin

#mstatspop v.0.1beta (20170224)

#zlib 1.2.8 installation (dependency)
#
#mkdir -p ./zlib
#wget http://zlib.net/zlib-1.3.tar.gz -P ./zlib
#tar -zxvf ./zlib/zlib-1.3.tar.gz -C ./zlib
#rm ./zlib/zlib-1.3.tar.gz
#cd ./zlib/zlib-1.3
#./configure
#make
#sudo make install

#gsl installation (dependency)
#mkdir -p /tmp/gsl && \
#    curl -o /tmp/gsl-2.2.tar.gz ftp://ftp.gnu.org/gnu/gsl/gsl-2.2.tar.gz -LOk && \
#    tar -zxvf /tmp/gsl-2.2.tar.gz -C /tmp/gsl && \
#    rm /tmp/gsl-2.2.tar.gz && \
#    cd /tmp/gsl/gsl-2.2 && \
#    ./configure && \
#    make && \
#    sudo make install

#perhaps is necessary to include a symbolic link for finding gsl libraries:
#for example:
#sudo ln -s /home/gonza057/anaconda3/lib/libgsl.a /usr/local/lib/libgsl.a

#To compile:
#gcc ./sources/*.c -lm -o ./bin/mstatspop -Wall -O3 -lz
#OR (IN CASE USING OPTIMAL TESTS include the GSL Scientific library, downloading from http://www.gnu.org/software/gsl/)
#for osx 10.12 add:-I/usr/local/include -L/usr/local/lib
#gcc ./sources/*.c -lgsl -lgslcblas -lm -o ./bin/mstatspop -Wall -DinGSL=1 -O3 -lz 
#or
#gcc ./sources/*.c -lgsl -lgslcblas -lm -o ./bin/mstatspop -Wall -DinGSL=1 -O3 -lz -L/opt/homebrew/Cellar/gsl/2.7.1/lib/ /opt/homebrew/Cellar/gsl/2.7.1/lib/libgsl.a -I/opt/homebrew/Cellar/gsl/2.7.1/include
#for linux:
#gcc ./sources/*.c -lgsl -lgslcblas -lm -o ./bin/mstatspop -Wall -DinGSL=1 -O3 -lz
