References:
If gcc is too old, follow this if you have CentOS:
https://ahelpme.com/linux/centos7/how-to-install-new-gcc-and-development-tools-under-centos-7/


To install Spheral and test installation:

cd
mkdir Spheral
mkdir github-Spheral
cd github-Spheral
git clone https://github.com/jmikeowen/spheral
cd spheral/src/
scl enable devtoolset-7 bash
./boot
mkdir BUILD
cd BUILD
../configure --prefix=/home/adestefa/Spheral --with-opt=3 --with-compilers=gnu --with-dbc=none 
make -j 20
cd /home/adestefa/github-Spheral/spheral/tests
/home/adestefa/Spheral/bin/ats -n 20 -e /home/adestefa/Spheral/bin/python integration.ats



To install VisIt:

cd
mkdir download-visit
cd download-visit/
wget http://portal.nersc.gov/project/visit/releases/2.13.3/visit2_13_3.linux-x86_64-rhel7.tar.gz
wget http://portal.nersc.gov/project/visit/releases/2.13.3/visit-install2_13_3
chmod 755 visit-install2_13_3
./visit-install2_13_3 2.13.3 linux-x86_64-rhel7 /home/adestefa/visit-2.13.3
cd
echo "export PATH=$PATH:/home/adestefa/visit-2.13.3/bin" >> .bashrc
