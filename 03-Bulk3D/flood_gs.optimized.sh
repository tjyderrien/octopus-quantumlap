#!/usr/bin/env bash
# Massive GS submitter for QuantumLaP (careful with resources!)
# Usage examples:
#   ./flood_gs.sh                          # normal
#   DRY_RUN=1 ./flood_gs.sh                # do not submit, just stage & print
#   PARALLEL_JOBS=4 ./flood_gs.sh          # try concurrent submits if GNU parallel is available
#   SKIP_EXISTING=1 ./flood_gs.sh          # skip if target folder already exists

set -Eeuo pipefail

# -------- knobs --------
PARALLEL_JOBS=${PARALLEL_JOBS:=-1}   # >1 uses GNU parallel if available
SKIP_EXISTING=${SKIP_EXISTING:=1}   # skip runs that already have an output folder
DRY_RUN=${DRY_RUN:=1}               # passed through to run_gs.sh

# -------- environment --------
: "${QuantumLaPsources:?Set QuantumLaPsources in your environment}"
: "${QuantumLaPruns:?Set QuantumLaPruns in your environment}"

# shellcheck source=/dev/null
source "$QuantumLaPsources/octopus_slabs/machines.env"
CheckComputerName
LoadEnvLibraries

echo "===== QuantumLaP Massive GS submitter ====="
echo "WARNING: ~19,200 core-hours per job (estimate). Handle with care."
echo "This wraps:  $QuantumLaPsources/03-Bulk3D/run_gs.sh"
[ "$DRY_RUN" -eq 1 ] && echo ">>> DRY-RUN mode is ON (no submissions) <<<"
[ "$SKIP_EXISTING" -eq 1 ] && echo ">>> Will skip targets that already exist <<<"

# -------- user-configurable sweeps --------
NumNodes=8
Mission="Si_GS_angle"

# Sample
SampleLayers_list=(1)     # e.g., 6 8 10 ...
SampleXY_list=(1)         # e.g., 8 10 12 ...

# k-points
kPoints_list=(40)         # e.g., 4 16 40

# Mesh (Å). If you want multiple, add more entries (0.23 0.20 ...)
MeshSizeX_SI_list=(0.23)
# Z mesh follows X mesh unless you explicitly create a list and use it.

# Physics knobs
InitialTemperature="0e0"
EnableSmearing=0
SmearingTe_list=(0)
ExtraStates=0

TheoryLevel="dft"         # or "independent_particles"
Functional="TB09"         # or "LDA"
AngleIncidence_list=(0)   # e.g., 0 0.25 0.5 ...

if [ "$EnableSmearing" -eq 0 ]; then
  SmearingTe_list=(0)
  ExtraStates=0
fi

# -------- assemble job matrix --------
build_cmds() {
  for iSmearingTe in "${SmearingTe_list[@]}"; do
    for iAngle in "${AngleIncidence_list[@]}"; do
      for ik in "${kPoints_list[@]}"; do
        for iL in "${SampleLayers_list[@]}"; do
          for iXY in "${SampleXY_list[@]}"; do
            for iMeshX in "${MeshSizeX_SI_list[@]}"; do
              iMeshZ="$iMeshX"

              # Predict target folder (must match run_gs.sh logic) for optional skipping:
              # We recreate the mesh in Bohr with 6 decimals so SKIP_EXISTING is accurate.
              meshX_AU=$(printf "%s\n" "$iMeshX / 0.52917721067" | bc -l)
              meshZ_AU=$(printf "%s\n" "$iMeshZ / 0.52917721067" | bc -l)
              meshX_AU_FMT=$(printf "%.6f" "$meshX_AU")
              meshZ_AU_FMT=$(printf "%.6f" "$meshZ_AU")
              if [ "$EnableSmearing" -eq 1 ]; then
                target="${QuantumLaPruns}/03-Bulk3D/Si-${iXY}x${iXY}x${iL}-k-${ik}x${ik}x1-Angle${iAngle}-Mesh${meshX_AU_FMT}-Mesh${meshZ_AU_FMT}-T${InitialTemperature}-Te${iSmearingTe}eV"
              else
                target="${QuantumLaPruns}/03-Bulk3D/Si-${iXY}x${iXY}x${iL}-k-${ik}x${ik}x1-Angle${iAngle}-Mesh${meshX_AU_FMT}-Mesh${meshZ_AU_FMT}-T${InitialTemperature}"
              fi

              if [ "$SKIP_EXISTING" -eq 1 ] && [ -d "${target}" ]; then
                echo "echo \"[SKIP] ${target} exists\"" 
                continue
              fi

              # Build the actual command (quoted, positional args preserved)
              cat <<EOF
echo "** SUBMIT ** Z=${iL} XY=${iXY} k=${ik} MeshX=${iMeshX}Å Angle=${iAngle} Smear=${iSmearingTe}eV";
DRY_RUN=${DRY_RUN} "$QuantumLaPsources/03-Bulk3D/run_gs.optimized.sh" \
  "${iL}" "${iXY}" "${iMeshX}" "${iMeshZ}" "${ik}" \
  "${InitialTemperature}" "${TheoryLevel}" "${Functional}" \
  "${Mission}" "${NumNodes}" "${EnableSmearing}" "${iSmearingTe}" \
  "${ExtraStates}" "${iAngle}"
EOF
            done
          done
        done
      done
    done
  done
}

# -------- execute job matrix --------
if command -v parallel >/dev/null 2>&1 && [ "$PARALLEL_JOBS" -gt 1 ]; then
  echo "Launching with GNU parallel (-j $PARALLEL_JOBS)"
  build_cmds | parallel -u -j "$PARALLEL_JOBS" bash -c "{}"
else
  [ "$PARALLEL_JOBS" -gt 1 ] && echo "[WARN] GNU parallel not found; running sequentially."
  while IFS= read -r line; do
    bash -c "$line"
  done < <(build_cmds)
fi

