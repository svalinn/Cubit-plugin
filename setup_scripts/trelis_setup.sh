#!/bin/bash

#Trelis Installation
INSTALL_ROOT=$HOME/cnerg
cd $INSTALL_ROUTE
mkdir trelis
cd trelis
echo Enter netID:
read username
curl -1 -v --disable-epsv --ftp-skip-pasv-ip -u $username@wisc.edu --ftp-ssl \
 --output trelis.deb \
 ftp://ftp.box.com/CNERG/INTERNAL/Resources/Trelis/Trelis-16.3.3-Lin64.deb > trelis.deb;
    #Currently uses Trelis 16.3, since SDK appears incompatible with 16.4
dpkg -i trelis.deb #This installs Trelis 16.3, files can be found in /opt

cd $INSTALL_ROOT
