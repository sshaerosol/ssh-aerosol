#!/bin/bash

export ssh_root=`pwd`
garbagedir=${ssh_root}/compilogs/
arch_path='myssh_arch'
compil_mode=PROF
tmplab=`date +"%s"`

if [ ! -d ${garbagedir} ]; then
    mkdir ${garbagedir}
fi

#clear previous compilation
rm -f `find -iname "*.mod"`
rm -f `find -iname "*.a"`
rm -f `find -iname "*.o"`
rm -rf Makefile.hdr
rm -rf mymodules.sh

# Deal with prompt 
while (($# > 0))
do
    case $1 in
        "-h"|"--h"|"--help"|"-help") 
            echo "build-ssh.sh - Compile SSH-aerosol with your defined architecture to be used with CHIMERE model"
            echo "build-chimere.sh [options]"
            echo "      [ -h | --h | -help | --help ] : this help message"
            echo "      [--arch arch] : to choose target architecture (file myssh_arch/myssh-arch must exist - default : last architecture used)"
            echo "      [--avail] : to know available target architectures"
            exit ;;
        "--arch") arch=$2 ; arch_defined=true ; shift ; shift ;;
        "--avail") ls ${arch_path}/mychimere-* | cut -d"-" -f2 | sed 's/^/    /'  ; exit ;;
        *) code=$1 ; shift ;;
    esac
done
compil=${arch##*.}
export my_mode=${compil_mode}
if test -f myssh_arch/makefiles.hdr/Makefile.hdr.${compil}-64-ompi; then
    ln -sfn myssh_arch/makefiles.hdr/Makefile.hdr.${compil}-64-ompi Makefile.hdr
else
    echo "Error in architecture definition. The file ./makefiles.hdr/Makefile.hdr.${compil}-64-ompi does not exist. Please check --arch [architecture choice]"
    exit 1
fi

if test -f myssh_arch/mychimere-${arch}; then
    ln -fsn myssh_arch/mychimere-${arch} mymodules.sh
else
    echo "Error in architecture definition. The file myssh_arch/mychimere-${arch} does not exist."
    exit 1
fi


echo "================================================="
echo "   Compiling SSH-aerosol..."
echo "================================================="
cd ${ssh_root}/src
ln -sfn ../mymodules.sh ./
ln -sfn ../Makefile.hdr Makefile.hdr
make clean > ${garbagedir}/make.ssh-aerosol.${tmplab}.log 2>&1
cd ../
./make.sh  > ${garbagedir}/make.ssh-aerosol.${tmplab}.log 2>&1
if [ $? -ne 0 ] ; then
    echo
    echo "================================================="
    tail -20 ${garbagedir}/make.ssh-aerosol.${tmplab}.log
    echo "================================================="
    echo -e "\033[01;31;1m SSH-aerosol compilation aborted \033[m"
    echo "Check file ${garbagedir}/make.ssh-aerosol.${tmplab}.log"
    echo
    exit 1
fi


N_warnings_gfortran=`grep -c "Warning" ${garbagedir}/make.ssh-aerosol.${tmplab}.log`
N_warnings_ifort=`grep -c "remark" ${garbagedir}/make.ssh-aerosol.${tmplab}.log`
N_warnings=$(($N_warnings_ifort+$N_warnings_gfortran))

echo -e " \033[1;42m  SSH-aerosol compilation OK \033[0m"

if [ ${N_warnings} -ne "0" ] ; then
    echo ${N_warnings}" Warnings in SSH-aerosol compilation. Check file ${garbagedir}/make.ssh-aerosol.${tmplab}.log"
fi

exit 0
