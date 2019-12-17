#!/bin/bash

echo "Count iNodes in subfolders. Useful for finding out who is taking all the space."
find | cut -d/ -f2 | uniq -c | sort -n
