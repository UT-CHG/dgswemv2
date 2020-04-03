#!/bin/bash

if [ "$#" -eq "0" ]; then
    cmd="/home/bremerm31/clang+llvm-6.0.1-x86_64-linux-gnu-ubuntu-16.04/bin/clang-format -i -style=file */*.cpp */*.hpp */*.tpp */*/*.cpp */*/*.hpp */*/*.tpp */*/*/*.cpp */*/*/*.hpp */*/*/*.tpp */*/*/*/*.cpp */*/*/*/*.hpp */*/*/*/*.tpp */*/*/*/*/*.cpp */*/*/*/*/*.hpp */*/*/*/*/*.tpp */*/*/*/*/*/*.cpp */*/*/*/*/*/*.hpp */*/*/*/*/*/*.tpp"
    echo ${cmd}
    ${cmd}
else
    for var in "$@"
    do
	cmd="clang-format-5.0 -i -style=file ${var}"
	echo ${cmd}
	${cmd}
    done
fi
