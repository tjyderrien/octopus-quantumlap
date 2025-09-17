#!/usr/bin/env bash

# Optimized by ChatGPT v5 on Sept 17th 2025. 
# Testing status: none. 

# Input (positional for backward compatibility):
#   $1: SampleLayers
#   $2: SampleXY
#   $3: MeshSizeXInput_SI (Å)
#   $4: MeshSizeZInput_SI (Å)
#   $5: kPoints
#   $6: InitialTemperature (string, e.g. 0e0)
#   $7: TheoryLevel ("dft" | "independent_particles")
#   $8: Functional ("LDA" | "TB09")
#   $9: Mission
#  $10: NumNodes
#  $11: EnableSmearing (0|1)
#  $12: SmearingTe (eV)
#  $13: ExtraStates (int)
#  $14: AngleIncidence (float)

set -Eeuo pipefail

# -------- configuration / defaults --------
DRY_RUN=${DRY_RUN:=-0}   # override with: DRY_RUN=1 ./run_gs.sh ...
BOHR_PER_ANG=0.52917721067

# -------- helpers --------
die() { printf "ERROR: %s\n" "$*" >&2; exit 1; }
info(){ printf "[INFO] %s\n" "$*"; }
warn(){ printf "[WARN] %s\n" "$*" >&2; }

calc() {
  # bc is POSIX-common and avoids starting Python processes repeatedly
  # usage: calc "3.0/2"
  printf "%s\n" "$1" | bc -l
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command '$1' not found in PATH."
}

require_file() {
  [ -f "$1" ] || die "Required file not found: $1"
}

# -------- environment --------
: "${QuantumLaPsources:?Set QuantumLaPsources in your environment}"
: "${QuantumLaPruns:?Set QuantumLaPruns in your environment}"

require_file "$QuantumLaPsources/octopus_slabs/machines.env"
# shellcheck source=/dev/null
source "$QuantumLaPsources/octopus_slabs/machines.env"

# These functions come from machines.env
type CheckComputerName >/dev/null 2>&1 || die "machines.env missing CheckComputerName()"
type LoadEnvLibraries >/dev/null 2>&1    || die "machines.env missing LoadEnvLibraries()"

CheckComputerName
LoadEnvLibraries

info "QuantumLaPsources: $QuantumLaPsources"
info "QuantumLaPruns:    $QuantumLaPruns"

# -------- parse positional inputs (kept for compatibility) --------
SampleLayers=${1:?Missing SampleLayers}
SampleXY=${2:?Missing SampleXY}
MeshSizeXInput_SI=${3:?Missing MeshSizeXInput_SI (Å)}
MeshSizeZInput_SI=${4:?Missing MeshSizeZInput_SI (Å)}
kPoints=${5:?Missing kPoints}
InitialTemperature=${6:?Missing InitialTemperature}
TheoryLevel=${7:?Missing TheoryLevel}
Functional=${8:?Missing Functional}
Mission=${9:?Missing Mission}
NumNodes=${10:?Missing NumNodes}
EnableSmearing=${11:?Missing EnableSmearing (0|1)}
SmearingTe=${12:?Missing SmearingTe}
ExtraStates=${13:?Missing ExtraStates}
AngleIncidence=${14:?Missing AngleIncidence}

# -------- validate & normalize --------
case "$TheoryLevel" in
  dft|independent_particles) ;;
  *) die "TheoryLevel must be 'dft' or 'independent_particles' (got: $TheoryLevel)";;
esac

case "$Functional" in
  LDA|TB09) ;;
  *) die "Functional must be 'LDA' or 'TB09' (got: $Functional)";;
# end_esac(){ :; } # shellcheck disable=SC2317
esac

case "$EnableSmearing" in
  0|1) ;;
  *) die "EnableSmearing must be 0 or 1 (got: $EnableSmearing)";;
esac

require_cmd bc
require_cmd sed
require_cmd awk

TheoryLevelPath=$([ "$TheoryLevel" = "independent_particles" ] && printf -- "-Indep_particles" || printf "")
FunctionalPath=$([ "$Functional" = "TB09" ] && printf -- "-TB09" || printf "")

# Convert mesh sizes from Å to Bohr (AU)
MeshSizeX_AU=$(calc "$MeshSizeXInput_SI / $BOHR_PER_ANG")
MeshSizeZ_AU=$(calc "$MeshSizeZInput_SI / $BOHR_PER_ANG")

# Truncate to a manageable number of decimals for path cleanliness
fmt_dec() { printf "%.6f" "$1"; }
MeshSizeX_AU_FMT=$(fmt_dec "$MeshSizeX_AU")
MeshSizeZ_AU_FMT=$(fmt_dec "$MeshSizeZ_AU")

# Build run path
if [ "$EnableSmearing" -eq 1 ]; then
  localPath="${QuantumLaPruns}/03-Bulk3D/Si-${SampleXY}x${SampleXY}x${SampleLayers}-k-${kPoints}x${kPoints}x1-Angle${AngleIncidence}-Mesh${MeshSizeX_AU_FMT}-Mesh${MeshSizeZ_AU_FMT}-T${InitialTemperature}-Te${SmearingTe}eV${TheoryLevelPath}${FunctionalPath}"
else
  localPath="${QuantumLaPruns}/03-Bulk3D/Si-${SampleXY}x${SampleXY}x${SampleLayers}-k-${kPoints}x${kPoints}x1-Angle${AngleIncidence}-Mesh${MeshSizeX_AU_FMT}-Mesh${MeshSizeZ_AU_FMT}-T${InitialTemperature}${TheoryLevelPath}${FunctionalPath}"
fi

info "Output folder: $localPath"
mkdir -p "$localPath"

# -------- input preparation --------
ReplaceInputFile() {
  info "Preparing 'inp' for ground state"
  [ -f inp ] && mv -f inp inp.bak

  require_file "$QuantumLaPsources/03-Bulk3D/inp-gs"
  cp "$QuantumLaPsources/03-Bulk3D/inp-gs" inp

  case "$Functional" in
    TB09)
      {
        echo "XCFunctional = mgga_x_tb09 + lda_c_pw"
        echo "FilterPotentials = filter_none"
      } >> inp
      ;;
    LDA)
      echo "XCFunctional = lda_x + lda_c_pz_mod" >> inp
      ;;
  esac
}

CustomizeInput() {
  # Replace placeholders in a single sed call where possible
  # Expecting these tokens to exist in the template: SizeOfTheBox, NumberOfBlocks, MeshSizeX_AU, MeshSizeZ_AU, kPointsValue, theta1_inp
  sed -i.bak \
    -e "s/SizeOfTheBox/${SampleXY}/g" \
    -e "s/NumberOfBlocks/${SampleLayers}/g" \
    -e "s/MeshSizeX_AU/${MeshSizeX_AU_FMT}/g" \
    -e "s/MeshSizeZ_AU/${MeshSizeZ_AU_FMT}/g" \
    -e "s/kPointsValue/${kPoints}/g" \
    -e "s/theta1_inp/${AngleIncidence}/g" \
    inp

  if [ "$TheoryLevel" = "independent_particles" ]; then
    warn "independent_particles selected. (Note: TDFreezeHxc is not enabled here.)"
    echo "TheoryLevel = independent_particles" >> inp
  fi

  if [ "$EnableSmearing" -eq 1 ]; then
    {
      echo "Smearing = ${SmearingTe} * eV"
      echo "SmearingFunction = fermi_dirac"
      echo "ExtraStates = ${ExtraStates}"
    } >> inp
  fi
}

SubmitTheJob() {
  info "Staging submission files"
  cp -u inp "$localPath/"
  cp "$QuantumLaPsources/03-Bulk3D/${SubmitScript}" "$localPath/"

  # Patch submission script
  USER_EMAIL="$(git config --get user.email || true)"
  [ -z "${USER_EMAIL}" ] && USER_EMAIL="noreply@example.com" && warn "git user.email not set; using ${USER_EMAIL}"

  sed -i.bak \
    -e "s/Mission/${Mission}/g" \
    -e "s/NumNodes/${NumNodes}/g" \
    -e "s/USER_EMAIL/${USER_EMAIL}/g" \
    "$localPath/${SubmitScript}"

  if [ "$DRY_RUN" -eq 1 ]; then
    warn "DRY-RUN enabled: not submitting."
    info  "${ComputerName}: ${SubmitCommand} ${SubmitScript}"
    return 0
  fi

  info "Submitting job via: ${SubmitCommand} ${SubmitScript}"
  (
    cd "$localPath"
    set +e
    "${SubmitCommand}" "${SubmitScript}" > submit.log 2>&1
    rc=$?
    set -e
    [ $rc -eq 0 ] || {
      warn "Submission command exited with code $rc"
      die "See $localPath/submit.log for details."
    }
    # Try to extract a job number (machine dependent)
    jobnumber="$(awk '{for(i=1;i<=NF;i++) if ($i ~ /^[0-9]+$/) print $i}' submit.log | tail -n1 || true)"
    if [ -n "$jobnumber" ]; then
      cp -f submit.log "${jobnumber}.log" || true
      info "Job number: $jobnumber"
    else
      warn "Could not detect job number from submit.log"
    fi
  )
}

# -------- run pipeline --------
info "Sample: Z=${SampleLayers} layers, XY=${SampleXY} atoms"
info "Mesh (Å): X=${MeshSizeXInput_SI}, Z=${MeshSizeZInput_SI}  ->  (Bohr): X=${MeshSizeX_AU_FMT}, Z=${MeshSizeZ_AU_FMT}"
ReplaceInputFile
CustomizeInput
SubmitTheJob

