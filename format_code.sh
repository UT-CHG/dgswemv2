#!/bin/bash

if [ "$#" -eq "0" ]; then
    cmd="clang-format -i -style=file */*.cpp */*.hpp */*.tpp */*/*.cpp */*/*.hpp */*/*.tpp */*/*/*.cpp */*/*/*.hpp */*/*/*.tpp */*/*/*/*.cpp */*/*/*/*.hpp */*/*/*/*.tpp */*/*/*/*/*.cpp */*/*/*/*/*.hpp */*/*/*/*/*.tpp */*/*/*/*/*/*.cpp */*/*/*/*/*/*.hpp */*/*/*/*/*/*.tpp"
    echo ${cmd}
    ${cmd}
else
    for var in "$@"
    do
	cmd="clang-format -i -style=file ${var}"
	echo ${cmd}
	${cmd}
    done
fi
