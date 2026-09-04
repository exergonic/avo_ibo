#!/bin/bash
# Regenerate all example outputs with the current pipeline at the
# project standard level (wB97X-D/def2-TZVP), charges per molecule.
cd "$(dirname "$0")/.." || exit 1
run_one() {  # $1 = xyz stem, $2 = charge
  echo "=== examples/$1.xyz (charge $2) $(date +%H:%M:%S)"
  env -u PYTHONPATH .pixi/envs/default/python.exe -m avogadro_ibo \
    --method wB97X-D --basis def2-TZVP --charge "$2" --spin 1 "examples/$1.xyz"
}
for m in allene ammonia benzene cyclooctatetraene cyclopropane diborane ethene malonaldehyde_enol methane methylamine ozone SO3 water zncl2; do
  run_one "$m" 0
done
for m in cyclopropenium ethylium tbutyl; do
  run_one "$m" 1
done
for m in cyclopropenyl_anion_planar cyclopropenyl_anion_nonplanar; do
  run_one "$m" -1
done
echo "ALL DONE $(date +%H:%M:%S)"
