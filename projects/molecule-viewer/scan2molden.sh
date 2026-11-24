#!/usr/bin/env bash
# scan2molden.sh - turn a Gaussian relaxed-scan .log into a file Molden can open,
# keeping only the converged geometry of every scan point (with its energy).
#
# Usage:
#   ./scan2molden.sh scan.log                  -> writes scan.molden
#   ./scan2molden.sh scan.log -x 8,12          -> leave out scan points 8 and 12 (outliers)
#   ./scan2molden.sh scan.log --xyz            -> writes scan-points.xyz (plain multi-frame .xyz) instead
#   ./scan2molden.sh scan.log -o other.molden  -> choose the output name
#
# Notes:
#   - A scan with N steps has N+1 points (the start structure plus N steps).
#   - Only the first job in the file is read; --Link1-- parts are ignored.
#   - Points where Gaussian printed "Optimization stopped" are kept but marked in the
#     comment line and listed as a warning, so you can decide to leave them out with -x.
#   - Needs only bash and awk (both are on every Linux cluster).

set -euo pipefail

usage() { sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

log=""; out=""; exclude=""; format="molden"
while [ $# -gt 0 ]; do
  case "$1" in
    -h|--help) usage 0 ;;
    -x|--exclude) exclude="${2:-}"; shift 2 ;;
    -o|--output) out="${2:-}"; shift 2 ;;
    --xyz) format="xyz"; shift ;;
    -*) echo "Unknown option: $1" >&2; usage 1 ;;
    *) if [ -z "$log" ]; then log="$1"; else echo "Only one log file at a time." >&2; exit 1; fi; shift ;;
  esac
done

[ -n "$log" ] || usage 1
[ -f "$log" ] || { echo "File not found: $log" >&2; exit 1; }
if [ -n "$exclude" ] && ! [[ "$exclude" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
  echo "Use -x with point numbers separated by commas, for example: -x 8,12" >&2; exit 1
fi

base="${log%.*}"
if [ -z "$out" ]; then
  if [ "$format" = "molden" ]; then out="$base.molden"; else out="$base-points.xyz"; fi
fi

# Same orientation everywhere: prefer Input orientation (what you gave Gaussian)
if grep -q "Input orientation:" "$log"; then orient="Input orientation:"
elif grep -q "Z-Matrix orientation:" "$log"; then orient="Z-Matrix orientation:"
elif grep -q "Standard orientation:" "$log"; then orient="Standard orientation:"
else echo "No geometries found in $log - is this a Gaussian output file?" >&2; exit 1
fi

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

awk -v orient="$orient" -v exclude="$exclude" -v format="$format" -v src="$(basename "$log")" -v warnfile="$tmp" '
BEGIN {
  split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn", sym, " ")
  n = split(exclude, ex, ","); for (i = 1; i <= n; i++) skip[ex[i] + 0] = 1
  nat = 0; ns = 0; np = 0; havegeom = 0
}
index($0, orient) {
  getline; getline; getline; getline            # dashes, two header lines, dashes
  nat = 0
  while ((getline line) > 0) {
    if (line ~ /^ *-----/) break
    c = split(line, p, " ")
    if (c < 5) continue
    nat++
    el = (p[2] + 0 >= 1 && p[2] + 0 <= 86) ? sym[p[2] + 0] : "X"
    geom[nat] = sprintf("%-2s %14.8f %14.8f %14.8f", el, p[c - 2], p[c - 1], p[c])
  }
  havegeom = 1
  next
}
/SCF Done:/ && havegeom {
  for (i = 1; i <= NF; i++) if ($i == "=") e = $(i + 1)
  ns++; stepE[ns] = e; stepN[ns] = nat
  for (k = 1; k <= nat; k++) stepG[ns, k] = geom[k]
  next
}
/Optimization completed/ { if (ns) { np++; pt[np] = ns; ptStat[np] = "converged" } ; next }
/Optimization stopped/   { if (ns) { np++; pt[np] = ns; ptStat[np] = "NOT converged" } ; next }
/Normal termination of Gaussian/ { exit }
END {
  if (np == 0) { print "NOPOINTS" > warnfile; exit }
  kept = 0
  for (j = 1; j <= np; j++) if (!(j in skip)) keep[++kept] = j
  if (kept == 0) { print "ALLSKIPPED" > warnfile; exit }
  if (format == "molden") {
    print "[Molden Format]"
    print "[N_GEO]"
    print kept
    print "[GEOCONV]"
    print "energy"
    for (q = 1; q <= kept; q++) print stepE[pt[keep[q]]]
    print "[GEOMETRIES] XYZ"
  }
  for (q = 1; q <= kept; q++) {
    j = keep[q]; s = pt[j]
    print stepN[s]
    printf "Energy: %s Eh | scan point %d of %d | %s | %s\n", stepE[s], j, np, ptStat[j], src
    for (k = 1; k <= stepN[s]; k++) print stepG[s, k]
  }
  print "POINTS " np " " kept > warnfile
  for (j = 1; j <= np; j++) if (ptStat[j] != "converged") print "STOPPED " j > warnfile
  for (j in skip) if (j + 0 > np) print "NOSUCH " j > warnfile
}
' "$log" > "$out"

if grep -q "^NOPOINTS" "$tmp"; then
  rm -f "$out"; echo "No finished scan points found in $log (no 'Optimization completed')." >&2; exit 1
fi
if grep -q "^ALLSKIPPED" "$tmp"; then
  rm -f "$out"; echo "All points were left out with -x, nothing to write." >&2; exit 1
fi

read -r _ total kept < <(grep "^POINTS" "$tmp")
echo "Read $total scan points from $log, wrote $kept to $out"
while read -r _ j; do echo "Warning: scan point $j did NOT converge (Optimization stopped). Leave it out with -x $j" >&2; done < <(grep "^STOPPED" "$tmp" || true)
while read -r _ j; do echo "Warning: there is no scan point $j (only $total points), -x $j was ignored." >&2; done < <(grep "^NOSUCH" "$tmp" || true)
