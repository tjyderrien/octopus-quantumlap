#!/bin/bash

python_calc () {
    python -c "print ($@)"
}

list_folders () {
    find . -maxdepth 4 -mindepth 0 -type d
}



