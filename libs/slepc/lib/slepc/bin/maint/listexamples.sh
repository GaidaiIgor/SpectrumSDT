#!/bin/bash
git ls-files -- src/*/tutorials/ex[1-9]* | awk -F/ '{printf "%-12s %s\n", $4, $2}' | sort -n -k 1.3
