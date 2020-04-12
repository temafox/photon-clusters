#!/bin/sh

if [[ -a src/algorithm/example/tmvaClassificationApplication.cpp ]]; then
	root -l src/algorithm/example/tmvaClassificationApplication.cpp
else
	echo 'Error: no classification application script found'
fi

