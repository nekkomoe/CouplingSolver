#!/bin/bash

# 会把$default_path复制到$cp2dirname
cp2dirname=$1"_"$2
# 如果cp2dirname不存在，或者该目录已存在，则退出
if [ -z $cp2dirname ] || [ -d $cp2dirname ]; then
    echo "Usage: $0 <cp2dirname> <index>"
    echo -e "dir \033[36m$cp2dirname\033[0m should not exist"
    exit 1
fi
# 如果$2不是>=1的整数，则退出
if ! [[ $2 =~ ^[1-9][0-9]*$ ]]; then
    echo "Usage: $0 <cp2dirname> <index>"
    echo -e "index \033[36m$2\033[0m should be a positive integer without leading zeros"
    exit 1
fi

# echo $flag
# exit 0

N=100
Nt=1000
dt=0.01
Tm=$(echo "$Nt * $dt" | bc)

Phi0=100.0
Ts0=1000
Tc0=290
U0=5.0

file_prepath=
# 如果$2是1，则要求$3是一个存在的文件
if [ $2 -eq 1 ]; then
    if [ -z $3 ]; then
        echo "Usage: $0 <cp2dirname> <index> <file_prepath>"
        echo -e "param \033[36m<file_prepath>\033[0m is empty"
        exit 1
    fi
    if [ ! -d $3 ]; then
        echo "Usage: $0 <cp2dirname> <index> <file_prepath>"
        echo -e "dir \033[36m$3\033[0m not found"
        exit 1
    fi
    file_prepath=$3
else
    file_prepath=$1"_"`expr $2 - 1`
fi

echo -e "load from \033[36m$file_prepath\033[0m"

file_pf=$file_prepath"/pf_final.bin"
file_pc=$file_prepath"/lc_final.bin"
if [ ! -f $file_pf ] || [ ! -f $file_pc ]; then
    echo "File not found: $file_pf or $file_pc"
    exit 1
fi

default_path="data"
default_exec="std"
mkdir -p $default_path
make clean
make \
  "CXXFLAG_DEF=$flag" $default_exec

# 如果编译失败，则退出
if [ ! -f $default_exec ]; then
    echo "Compile failed"
    exit 1
fi

mpirun -n $pnum ./std \
    -file_pf $file_pf \
    -file_lc $file_pc \
    -Nr $N -Nz $N   \
    -Nt $Nt -Tm $Tm \
    -Tc0 $Tc0       \
    -U0 $U0         \
    -ksp_type gmres \
    -pc_type bjacobi $para

cp -r data $cp2dirname
echo -e "save to \033[36m$cp2dirname\033[0m"
