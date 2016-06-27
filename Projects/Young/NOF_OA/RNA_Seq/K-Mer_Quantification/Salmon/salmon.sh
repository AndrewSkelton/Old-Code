#!/bin/bash
#--extraSensitive
salmon quant \
				    -g ${5} \
			      -p 2 \
			      -i ${3} \
			      -l 'IU' \
			      -1 <(gunzip -c ${1}) \
			      -2 <(gunzip -c ${2}) \
			      -o ${4} &
