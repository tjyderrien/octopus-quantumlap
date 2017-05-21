#!/bin/bash

python_calc () {
    python -c "print ($@)"
}

list_folders () {
    find . -type d
}

