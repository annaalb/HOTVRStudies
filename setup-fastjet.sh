cd /nfs/dust/cms/user/albrecha/uhh2_102X_v2/fastjet-contrib/
make -j9
make check
make install
make fragile-shared
make fragile-shared-install
cd /nfs/dust/cms/user/albrecha/uhh2_102X_v2/CMSSW_10_2_17/src/UHH2/HOTVRStudies
make -j9
