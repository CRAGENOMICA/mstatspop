#mstatspop v.0.1beta (20170224)

#zlib 1.2.8 installation (dependency)
#
#mkdir -p ./zlib
#wget http://zlib.net/zlib-1.2.8.tar.gz -P ./zlib
#tar -zxvf ./zlib/zlib-1.2.8.tar.gz -C ./zlib
#rm ./zlib/zlib-1.2.8.tar.gz
#cd ./zlib/zlib-1.2.8
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

#To compile:
#gcc ./sources/*.c -lm -o ./bin/mstatspop -Wall -O3 -lz
#OR (IN CASE USING OPTIMAL TESTS include the GSL Scientific library, downloading from http://www.gnu.org/software/gsl/)
#for osx 10.12 add:-I/usr/local/include -L/usr/local/lib
gcc ./sources/*.c -lgsl -lgslcblas -lm -o ./bin/mstatspop -Wall -DinGSL=1 -O3 -L/usr/local/lib -I/usr/local/include /usr/local/lib/libgsl.a -lz
#for linux:
#gcc ./sources/*.c -lgsl -lgslcblas -lm -o ./bin/mstatspop -Wall -DinGSL=1 -O3 -lz
