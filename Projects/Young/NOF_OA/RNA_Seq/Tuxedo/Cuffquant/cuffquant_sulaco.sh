#!/bin/bash
cuffquant -o ${1} \
          -p 2 \
          -b ${4} \
          --max-bundle-frags 500000000 \
          ${5} \
          ${6}

mv ${1}/abundances.cxb ${2}/abundances/${3}.cbx

