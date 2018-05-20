#!/bin/bash

python_calc () {
    python -c "print ($@)"
}

list_folders () {
    find . -maxdepth 3 -mindepth 0 -type d
}

