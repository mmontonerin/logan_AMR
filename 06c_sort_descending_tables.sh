#!/bin/bash

for table in ./data/organism_counts*full*.csv
do
    filename=$(basename ${table} ".csv")
    sort -t, -k2,2nr ${table} > ./data/sorted_${filename}.csv
done