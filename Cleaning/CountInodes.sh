#!/bin/bash
find | cut -d/ -f2 | uniq -c | sort -n
