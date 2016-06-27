#!/bin/bash

kallisto quant -i ${1} -o ${4} ${2} ${3}

mv ${4}/abundance.txt ${5}/${6}.txt
