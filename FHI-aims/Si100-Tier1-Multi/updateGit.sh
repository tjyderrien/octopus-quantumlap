#!/bin/bash
echo "[Git] Adding band structures to repository..."
git add `find | grep 'band'`

echo "[Git] Adding DOS to repository..."
git add `find | grep 'KS_DOS'`

echo "[Git] Adding geometries to repository..."
git add `find | grep 'geometry.in'`

echo "[Git] Adding controls of FHI-aims to repository..."
git add `find | grep 'control.in'`

echo "[Git] Adding FHI-aims simulation results to repository..."
git add `find | grep 'aims.out'`

echo "[Git] Adding pictures generated from plotting DOS..."
git add `find | grep 'DOS.eps'`
git add `find | grep 'DOS.png'`

echo "git commit -m 'Put here your message for committing"
echo "git push"
