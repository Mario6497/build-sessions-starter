#!/usr/bin/env python3
"""molstudio.py - Molecule studio, a local Molden-like viewer (build step 2).

Opens Gaussian .log/.out, ORCA .out, ORCA _trj.xyz and .xyz files (single or
multi-frame). Relaxed scans show only the converged point of every scan step;
optimisations and TS searches show every step. Energy plot on the left,
rotatable 3D structure on the right.
Vibrations tab (Gaussian / ORCA frequency jobs): frequency list, animated
normal modes with arrows, IR spectrum, and forward/backward structures
displaced along a mode (e.g. the imaginary mode of a TS) saved as .xyz.

Start (Ubuntu):
    python3 molstudio.py              opens the made-up sample scan
    python3 molstudio.py scan.log     opens your file

Keys: Ctrl+O open · Ctrl+Shift+S save this structure · Ctrl+S save all points ·Left/Right step through points (or modes)
      Delete remove marked points · Ctrl+Z undo
      1 ball and stick · 2 sticks · 3 spacefill · 4 wireframe · N atom numbers on/off
Needs only the Python standard library with tkinter (sudo apt install python3-tk).
"""
import glob
import math
import os
import re
import sys
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

APP = "Molecule studio"

# ---------------------------------------------------------------- element data
SYMBOLS = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
           "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr "
           "Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn").split()
COV = {"H": .31, "He": .28, "Li": 1.28, "Be": .96, "B": .84, "C": .76, "N": .71, "O": .66, "F": .57,
       "Ne": .58, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
       "Ar": 1.06, "K": 2.03, "Ca": 1.76, "Ti": 1.6, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32,
       "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.2, "As": 1.19, "Se": 1.2,
       "Br": 1.2, "Kr": 1.16, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44, "Sn": 1.39,
       "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.4, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Pb": 1.46,
       "Bi": 1.48}
VDW = {"H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52,
       "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80,
       "Cl": 1.75, "Ar": 1.88, "K": 2.75, "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.85,
       "Kr": 2.02, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16}
COLORS = {"H": "#ffffff", "C": "#909090", "N": "#3050f8", "O": "#ff0d0d", "F": "#90e050",
          "Cl": "#1ff01f", "Br": "#a62929", "I": "#940094", "S": "#ffff30", "P": "#ff8000",
          "B": "#ffb5b5", "Si": "#f0c8a0", "Fe": "#e06633", "Cu": "#c88033", "Zn": "#7d80b0",
          "Pd": "#006985", "Pt": "#d0d0e0", "Au": "#ffd123", "Na": "#ab5cf2", "K": "#8f40d4",
          "Li": "#cc80ff", "Mg": "#8aff00"}
HARTREE_TO_KCAL = 627.5095
FLOAT = r"-?\d+\.\d+(?:[eEdD][-+]?\d+)?"

MODELS = (("ball", "Ball and stick"), ("stick", "Sticks"), ("space", "Spacefill"), ("wire", "Wireframe"))
STICK_R = 0.16       # stick radius in Å


def symbol_of(z):
    return SYMBOLS[z] if 0 < z < len(SYMBOLS) else "X"


def norm_symbol(s):
    if s.isdigit():
        return symbol_of(int(s))
    m = re.match(r"[A-Za-z]{1,2}", s)
    if not m:
        return "X"
    t = m.group(0)[0].upper() + m.group(0)[1:].lower()
    return t if t in SYMBOLS else t[0]


def color_of(el):
    return COLORS.get(el, "#ff1493")


def atom_radius(el, model):
    if model == "ball":
        return 0.2 + 0.25 * COV.get(el, 1.4)
    if model == "stick":
        return STICK_R
    if model == "space":
        return VDW.get(el, 2.0)
    return 0.0


class Step:
    def __init__(self, atoms, energy, status=""):
        self.atoms = atoms          # list of (element, x, y, z)
        self.energy = energy        # Hartree or None
        self.status = status        # "", "converged", "not converged"


class Result:
    """What a parser found in a file."""

    def __init__(self, kind, frames, scan_coord=None, program=""):
        self.kind = kind            # scan, opt, steps, xyz, single
        self.frames = frames
        self.scan_coord = scan_coord  # (type letter, [0-based atom indices]) or None
        self.program = program      # gaussian, orca, xyz
        self.vib = None             # Vibrations or None


class Mode:
    def __init__(self, label, freq, ir, vec):
        self.label = label          # mode number as the program prints it
        self.freq = freq            # cm-1, negative = imaginary
        self.ir = ir                # km/mol
        self.vec = vec              # list of (dx, dy, dz) per atom


class Vibrations:
    def __init__(self, atoms, modes):
        self.atoms = atoms          # geometry the modes belong to
        self.modes = modes


def read_orient(lines, i):
    """Gaussian orientation table starting at line i -> list of atoms."""
    atoms = []
    j = i + 5
    while j < len(lines) and not re.match(r"^\s*-{5,}", lines[j]):
        p = lines[j].split()
        if len(p) >= 5:
            atoms.append((symbol_of(int(p[1])), float(p[-3]), float(p[-2]), float(p[-1])))
        j += 1
    return atoms


def parse_gaussian_freq(text):
    """Last frequency section of a Gaussian output, or None."""
    lines = text.splitlines()
    start = None
    for k, l in enumerate(lines):
        if "Harmonic frequencies" in l:
            start = k
    if start is None:
        return None
    geom = None
    for orient in ("Standard orientation:", "Input orientation:", "Z-Matrix orientation:"):
        idx = [k for k in range(start) if orient in lines[k]]
        if idx:
            geom = read_orient(lines, idx[-1])
            break
    if not geom:
        return None
    n = len(geom)
    modes = []
    i = start
    while i < len(lines):
        l = lines[i]
        if "- Thermochemistry -" in l or "Normal termination" in l:
            break
        if re.match(r"^\s*Frequencies --\s", l):
            freqs = [float(x) for x in l.split("--", 1)[1].split()]
            nums = lines[i - 2].split()
            labels = nums if len(nums) == len(freqs) and all(x.isdigit() for x in nums) \
                else [str(len(modes) + k + 1) for k in range(len(freqs))]
            ir = [0.0] * len(freqs)
            j = i + 1
            while j < len(lines) and not re.match(r"^\s*Atom\s+AN\s+X", lines[j]):
                if lines[j].strip().startswith("IR Inten"):
                    ir = [float(x) for x in lines[j].split("--", 1)[1].split()]
                j += 1
            vecs = [[] for _ in freqs]
            for a in range(n):
                p = lines[j + 1 + a].split()[2:]
                for m in range(len(freqs)):
                    vecs[m].append(tuple(float(x) for x in p[3 * m:3 * m + 3]))
            for m in range(len(freqs)):
                modes.append(Mode(labels[m], freqs[m], ir[m] if m < len(ir) else 0.0, vecs[m]))
            i = j + 1 + n
            continue
        i += 1
    return Vibrations(geom, modes) if modes else None


def parse_orca_freq(text):
    """Last frequency section of an ORCA output, or None."""
    lines = text.splitlines()
    vf = None
    for k, l in enumerate(lines):
        if "VIBRATIONAL FREQUENCIES" in l:
            vf = k
    if vf is None:
        return None
    geom = None
    for k in range(vf, -1, -1):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in lines[k]:
            geom = []
            for l in lines[k + 2:]:
                p = l.split()
                if len(p) < 4:
                    break
                geom.append((norm_symbol(p[0]), float(p[1]), float(p[2]), float(p[3])))
            break
    if not geom:
        return None
    n = len(geom)
    freqs, nm, irs = {}, None, {}
    for k in range(vf + 1, len(lines)):
        m = re.match(r"^\s*(\d+):\s+(-?\d+\.\d+)\s+cm\*\*-1", lines[k])
        if m:
            freqs[int(m.group(1))] = float(m.group(2))
        elif "NORMAL MODES" in lines[k]:
            nm = k
            break
    if nm is None:
        return None
    cols, cur = {}, None
    k = nm + 1
    while k < len(lines):
        l = lines[k]
        if "IR SPECTRUM" in l or "THERMOCHEMISTRY" in l:
            break
        t = l.split()
        if t and all(x.isdigit() for x in t):
            cur = [int(x) for x in t]
        elif cur and t and t[0].isdigit() and len(t) == len(cur) + 1:
            row = int(t[0])
            for c, v in zip(cur, t[1:]):
                cols.setdefault(c, {})[row] = float(v)
        k += 1
    if k < len(lines) and "IR SPECTRUM" in lines[k]:
        for l in lines[k + 1:]:
            m = re.match(r"^\s*(\d+):\s+(-?\d+\.\d+)\s+\S+\s+(-?\d+\.\d+)", l)
            if m:
                irs[int(m.group(1))] = float(m.group(3))
            elif irs and not l.strip():
                break
    modes = []
    for c in sorted(freqs):
        if abs(freqs[c]) < 1e-3 or c not in cols or len(cols[c]) < 3 * n:
            continue
        v = cols[c]
        modes.append(Mode(str(c), freqs[c], irs.get(c, 0.0),
                          [(v[3 * a], v[3 * a + 1], v[3 * a + 2]) for a in range(n)]))
    return Vibrations(geom, modes) if modes else None


def max_move(vec):
    return max((length((0,) + tuple(v)) for v in vec), default=0.0)


def displaced(atoms, vec, step):
    """Move along a mode so that the atom that moves most moves `step` Å (negative = backward)."""
    m = max_move(vec) or 1.0
    return [(el, x + step * v[0] / m, y + step * v[1] / m, z + step * v[2] / m)
            for (el, x, y, z), v in zip(atoms, vec)]


def ir_curve(modes, lo, hi, fwhm, npts):
    """Sum of Lorentz peaks; every peak is as high as its IR intensity."""
    g = (fwhm / 2.0) ** 2
    out = []
    for k in range(npts + 1):
        nu = lo + (hi - lo) * k / float(npts)
        out.append((nu, sum(m.ir * g / ((nu - m.freq) ** 2 + g) for m in modes)))
    return out


def pick_frames(steps, is_opt, is_scan, what):
    """Scans: only finished scan points. Optimisations: every step up to the end."""
    if is_scan:
        frames = [s for s in steps if s.status]
        if not frames:
            raise ValueError("This %s scan has no finished scan points yet." % what)
        return "scan", frames
    if is_opt:
        end = next((k for k, s in enumerate(steps) if s.status), None)
        return "opt", (steps[:end + 1] if end is not None else steps)
    if len(steps) > 1:
        return "steps", steps
    return "single", steps[-1:]


# ---------------------------------------------------------------- Gaussian
def parse_gaussian(text):
    lines = text.splitlines()
    # only the first job (ignore --Link1-- parts such as freq or a single point)
    for k, l in enumerate(lines):
        if "Normal termination of Gaussian" in l:
            lines = lines[:k + 1]
            break
    first_job = "\n".join(lines)
    orient = next((t for t in ("Input orientation:", "Z-Matrix orientation:", "Standard orientation:")
                   if t in first_job), None)
    if orient is None:
        raise ValueError("No geometry found in this Gaussian file.")

    route = ""
    for i, l in enumerate(lines):
        if re.match(r"^\s#", l):
            j = i
            while j < len(lines) and not re.match(r"^\s-{5,}", lines[j]):
                route += lines[j].strip()
                j += 1
            break
    is_opt = re.search(r"(^|[\s#/])opt\b", route, re.I) is not None

    scan_coord = None
    for i, l in enumerate(lines):
        if "ModRedundant input section has been read" in l:
            for j in range(i + 1, min(len(lines), i + 50)):
                m = re.match(r"^\s*([BADL])((?:\s+\d+){2,4})\s+S\s+\d+", lines[j], re.I)
                if m:
                    scan_coord = (m.group(1).upper(), [int(n) - 1 for n in m.group(2).split()])
                    break
                if not lines[j].strip():
                    break
            break
    is_scan = is_opt and (scan_coord is not None or "Potential Surface Scan" in first_job
                          or "scan point" in first_job)

    steps = []
    geom = None
    scf_re = re.compile(r"SCF Done:\s+E\([^)]*\)\s*=\s*(" + FLOAT + ")")
    i = 0
    while i < len(lines):
        l = lines[i]
        if orient in l:
            atoms = []
            j = i + 5
            while j < len(lines) and not re.match(r"^\s*-{5,}", lines[j]):
                p = lines[j].split()
                if len(p) >= 5:
                    atoms.append((symbol_of(int(p[1])), float(p[-3]), float(p[-2]), float(p[-1])))
                j += 1
            geom = atoms
            i = j + 1
            continue
        m = scf_re.search(l)
        if m and geom:
            steps.append(Step(geom, float(re.sub("[dD]", "e", m.group(1)))))
        elif "Optimization completed" in l and steps:
            steps[-1].status = "converged"
        elif "Optimization stopped" in l and steps:
            steps[-1].status = "not converged"
        i += 1
    if not steps:
        raise ValueError("Found geometries but no 'SCF Done' energies in this Gaussian file.")
    kind, frames = pick_frames(steps, is_opt, is_scan, "Gaussian")
    res = Result(kind, frames, scan_coord if kind == "scan" else None, "gaussian")
    res.vib = parse_gaussian_freq(text)      # whole file: freq often runs in a --Link1-- part
    return res


# ---------------------------------------------------------------- ORCA
def parse_orca(text):
    lines = text.splitlines()
    inp = []
    for l in lines:
        m = re.match(r"^\|\s*\d+>\s?(.*)$", l)
        if m:
            inp.append(m.group(1))
        elif inp and "END OF INPUT" in l:
            break
    inp_text = " ".join(inp)
    keywords = " ".join(x for x in inp if x.strip().startswith("!")).lower()
    is_opt = (re.search(r"\b(opt|optts|copt|zopt|tightopt|looseopt|verytightopt)\b", keywords) is not None
              or "Geometry Optimization Run" in text)
    scan_coord = None
    m = re.search(r"\bscan\b.*?\b([BAD])\s+(\d+)\s+(\d+)(?:\s+(\d+))?(?:\s+(\d+))?\s*=", inp_text, re.I)
    if m:
        t = m.group(1).upper()
        idx = [int(g) for g in m.groups()[1:] if g is not None][:{"B": 2, "A": 3, "D": 4}[t]]
        scan_coord = (t, idx)
    is_scan = is_opt and (scan_coord is not None or "RELAXED SURFACE SCAN STEP" in text)

    steps, geom, replace_next = [], None, False
    for i, l in enumerate(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in l:
            atoms = []
            j = i + 2
            while j < len(lines) and lines[j].strip():
                p = lines[j].split()
                try:
                    atoms.append((norm_symbol(p[0]), float(p[1]), float(p[2]), float(p[3])))
                except (ValueError, IndexError):
                    break
                j += 1
            if atoms:
                geom = atoms
        elif l.strip().startswith("FINAL SINGLE POINT ENERGY") and geom:
            try:
                e = float(l.split()[-1])
            except ValueError:
                continue
            st = Step(geom, e)
            if replace_next and steps:       # final evaluation of the converged structure
                st.status = steps[-1].status
                steps[-1] = st
                replace_next = False
            else:
                steps.append(st)
        elif "THE OPTIMIZATION HAS CONVERGED" in l and steps:
            steps[-1].status = "converged"
        elif "optimization did not converge" in l.lower() and steps:
            steps[-1].status = "not converged"
        elif "FINAL ENERGY EVALUATION AT THE STATIONARY POINT" in l:
            replace_next = True
    if not steps:
        raise ValueError("Found no 'FINAL SINGLE POINT ENERGY' with a geometry in this ORCA file.")
    kind, frames = pick_frames(steps, is_opt, is_scan, "ORCA")
    res = Result(kind, frames, scan_coord if kind == "scan" else None, "orca")
    res.vib = parse_orca_freq(text)
    return res


# ---------------------------------------------------------------- xyz
def parse_xyz(text):
    lines = text.splitlines()
    frames = []
    i = 0
    e_re1 = re.compile(r"(?:energy|\bE)\s*[:=]?\s*(" + FLOAT + ")", re.I)
    e_re2 = re.compile("(" + FLOAT + ")")
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        head = lines[i].strip()
        if not head.isdigit():
            raise ValueError("Line %d should be the number of atoms." % (i + 1))
        n = int(head)
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        atoms = []
        for k in range(n):
            p = lines[i + 2 + k].split() if i + 2 + k < len(lines) else []
            if len(p) < 4:
                raise ValueError("Frame %d has fewer atom lines than it says (%d)." % (len(frames) + 1, n))
            atoms.append((norm_symbol(p[0]), float(p[1]), float(p[2]), float(p[3])))
        m = e_re1.search(comment) or e_re2.search(comment)
        frames.append(Step(atoms, float(m.group(1)) if m else None))
        i += 2 + n
    if not frames:
        raise ValueError("No structures found in this .xyz file.")
    return Result("xyz" if len(frames) > 1 else "single", frames, None, "xyz")


def parse_any(name, text):
    if re.search(r"Entering Gaussian System|Gaussian, Inc\.|Input orientation:|Standard orientation:", text):
        return parse_gaussian(text)
    if re.search(r"\*\s+O\s+R\s+C\s+A\s+\*|ORCA TERMINATED|CARTESIAN COORDINATES \(ANGSTROEM\)", text):
        return parse_orca(text)
    if re.search(r"\.(mol|sdf|pdb)$", name, re.I):
        raise ValueError(".mol / .sdf / .pdb files are not part of this version.")
    first = next((l for l in text.splitlines() if l.strip()), "")
    if first.strip().isdigit():
        res = parse_xyz(text)
        if name.lower().endswith("_trj.xyz"):     # ORCA optimisation trajectory
            res.program = "orca"
            if res.kind == "xyz":
                res.kind = "opt"
        return res
    raise ValueError("Could not recognise this file. Supported: Gaussian .log/.out, ORCA .out, .xyz.")


# ---------------------------------------------------------------- geometry
def sub(a, b):
    return (a[1] - b[1], a[2] - b[2], a[3] - b[3])


def dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def length(u):
    return math.sqrt(dot(u, u))


def measure_value(atoms, idx):
    if any(i < 0 or i >= len(atoms) for i in idx):
        return None
    p = [atoms[i] for i in idx]
    if len(idx) == 2:
        return "Bond", length(sub(p[0], p[1])), "Å"
    if len(idx) == 3:
        u, v = sub(p[0], p[1]), sub(p[2], p[1])
        c = max(-1.0, min(1.0, dot(u, v) / (length(u) * length(v))))
        return "Angle", math.degrees(math.acos(c)), "°"
    if len(idx) == 4:
        b1, b2, b3 = sub(p[1], p[0]), sub(p[2], p[1]), sub(p[3], p[2])
        n1, n2 = cross(b1, b2), cross(b2, b3)
        b2hat = tuple(x / length(b2) for x in b2)
        return "Dihedral", math.degrees(math.atan2(dot(cross(n1, n2), b2hat), dot(n1, n2))), "°"
    return None


def bonds_of(atoms):
    out = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            lim = COV.get(atoms[i][0], 1.4) + COV.get(atoms[j][0], 1.4) + 0.4
            d = length(sub(atoms[i], atoms[j]))
            if 0.1 < d < lim:
                out.append((i, j))
    return out


def centroid(atoms):
    n = float(len(atoms) or 1)
    return (sum(a[1] for a in atoms) / n, sum(a[2] for a in atoms) / n, sum(a[3] for a in atoms) / n)


def mat_mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def rot_x(a):
    return [[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]]


def rot_y(a):
    return [[math.cos(a), 0, math.sin(a)], [0, 1, 0], [-math.sin(a), 0, math.cos(a)]]


def shade(hexcol, f):
    r, g, b = int(hexcol[1:3], 16), int(hexcol[3:5], 16), int(hexcol[5:7], 16)
    if f > 0:
        r, g, b = r + (255 - r) * f, g + (255 - g) * f, b + (255 - b) * f
    else:
        r, g, b = r * (1 + f), g * (1 + f), b * (1 + f)
    return "#%02x%02x%02x" % (int(r), int(g), int(b))


def nice_ticks(lo, hi, count):
    span = (hi - lo) or 1.0
    raw = span / count
    mag = 10 ** math.floor(math.log10(raw))
    r = raw / mag
    step = (10 if r >= 7.5 else 5 if r >= 3.5 else 2 if r >= 1.5 else 1) * mag
    t, v = [], math.ceil(lo / step - 1e-9) * step
    while v <= hi + step * 1e-6:
        t.append(round(v, 10))
        v += step
    return t


# ---------------------------------------------------------------- the app
class App:
    LABELS = {"scan": "Scan point", "opt": "Optimisation step", "xyz": "Frame",
              "steps": "Structure", "single": "Structure"}

    def __init__(self, root):
        self.root = root
        root.title(APP)
        root.geometry("1250x800")
        self.here = os.path.dirname(os.path.abspath(__file__))
        self.path, self.name, self.kind, self.frames, self.scan_coord = "", "", "", [], None
        self.program = ""
        self.removed, self.undo_stack, self.selected, self.current = set(), [], set(), 0
        self.rot, self.zoom, self.extent = self.default_rot(), 1.0, 5.0
        self.measure, self.bond_cache, self.projected, self.plot_pts = [], {}, [], []
        self.drag = None
        self.tree_busy = False
        self.redraw_pending = False
        self.show_labels = tk.BooleanVar(value=True)
        self.zero_based = tk.BooleanVar(value=False)
        self.model = tk.StringVar(value="ball")
        self.vib, self.vib_sel, self.vib_busy, self.spec_pts = None, None, False, []
        self.anim_on, self.anim_job, self.phase = False, None, 0.0
        self.amp = tk.DoubleVar(value=0.4)
        self.show_arrows = tk.BooleanVar(value=True)
        self.fwhm = tk.StringVar(value="20")
        self.disp_step = tk.StringVar(value="0.10")
        self.build_menu()
        self.build_ui()

    # ---------- menu
    def build_menu(self):
        bar = tk.Menu(self.root)
        fm = tk.Menu(bar, tearoff=0)
        fm.add_command(label="Open…", accelerator="Ctrl+O", command=self.open_dialog)
        sm = tk.Menu(fm, tearoff=0)
        for p in sorted(glob.glob(os.path.join(self.here, "sample-*"))):
            sm.add_command(label=os.path.basename(p), command=lambda p=p: self.load_path(p))
        fm.add_cascade(label="Open sample", menu=sm)
        fm.add_separator()
        fm.add_command(label="Save this structure as .xyz…", accelerator="Ctrl+Shift+S",
                       command=self.save_current)
        fm.add_command(label="Save all remaining points as .xyz…", accelerator="Ctrl+S", command=self.save_xyz)
        fm.add_separator()
        fm.add_command(label="Quit", accelerator="Ctrl+Q", command=self.root.destroy)
        bar.add_cascade(label="File", menu=fm)
        vm = tk.Menu(bar, tearoff=0)
        for k, (key, label) in enumerate(MODELS):
            vm.add_radiobutton(label=label, value=key, variable=self.model, accelerator=str(k + 1),
                               command=self.schedule_view)
        vm.add_separator()
        vm.add_checkbutton(label="Atom numbers", variable=self.show_labels, accelerator="N",
                           command=self.schedule_view)
        vm.add_checkbutton(label="Count atoms from 0 (ORCA style)", variable=self.zero_based,
                           command=self.numbering_changed)
        vm.add_separator()
        vm.add_command(label="Reset view", command=self.reset_view)
        vm.add_command(label="Clear measurement", command=self.clear_measure)
        bar.add_cascade(label="View", menu=vm)
        self.root.config(menu=bar)

    # ---------- UI
    def build_ui(self):
        top = ttk.Frame(self.root, padding=6)
        top.pack(fill="x")
        ttk.Button(top, text="Open… (Ctrl+O)", command=self.open_dialog).pack(side="left")
        self.b_remove = ttk.Button(top, text="Remove point", command=self.remove_action)
        self.b_remove.pack(side="left", padx=(12, 0))
        self.b_undo = ttk.Button(top, text="Undo", command=self.undo)
        self.b_undo.pack(side="left", padx=4)
        self.b_restore = ttk.Button(top, text="Bring all back", command=self.restore_all)
        self.b_restore.pack(side="left")
        ttk.Button(top, text="Save this structure (Ctrl+Shift+S)",
                   command=self.save_current).pack(side="left", padx=(12, 4))
        ttk.Button(top, text="Save all points (Ctrl+S)", command=self.save_xyz).pack(side="left")
        self.info = ttk.Label(top, text="", foreground="#555")
        self.info.pack(side="left", padx=12)

        paned = ttk.PanedWindow(self.root, orient="horizontal")
        paned.pack(fill="both", expand=True, padx=6, pady=(0, 6))
        left, right = ttk.Frame(paned), ttk.Frame(paned)
        paned.add(left, weight=1)
        paned.add(right, weight=1)

        self.nb = ttk.Notebook(left)
        self.nb.pack(fill="both", expand=True)
        ef, vf = ttk.Frame(self.nb, padding=4), ttk.Frame(self.nb, padding=4)
        self.energy_tab, self.vib_tab = ef, vf
        self.nb.add(ef, text="Energy")
        self.nb.add(vf, text="Vibrations")
        self.nb.bind("<<NotebookTabChanged>>", self.on_tab)
        self.build_vib_tab(vf)

        ttk.Label(ef, text="Click a point to see its structure · Ctrl+click marks several · "
                           "←/→ step · Delete removes · double-click a removed row to bring it back",
                  foreground="#666", wraplength=560).pack(anchor="w")
        self.plot = tk.Canvas(ef, height=320, bg="white", highlightthickness=0)
        self.plot.pack(fill="x", pady=4)
        self.plot.bind("<Button-1>", lambda e: self.plot_click(e, False))
        self.plot.bind("<Control-Button-1>", lambda e: self.plot_click(e, True))
        self.plot.bind("<Configure>", lambda e: self.draw_plot())

        tf = ttk.Frame(ef)
        tf.pack(fill="both", expand=True)
        cols = ("pt", "coord", "e", "de", "status")
        self.tree = ttk.Treeview(tf, columns=cols, show="headings", selectmode="extended")
        for c, t, w in zip(cols, ("Point", "Coordinate", "E (Eh)", "ΔE (kcal/mol)", "Status"),
                           (60, 110, 140, 110, 110)):
            self.tree.heading(c, text=t)
            self.tree.column(c, width=w, anchor="e")
        self.tree.tag_configure("removed", foreground="#a0a7b4")
        sb = ttk.Scrollbar(tf, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        self.tree.bind("<<TreeviewSelect>>", self.tree_select)
        self.tree.bind("<Double-1>", self.tree_double)

        self.title = ttk.Label(right, text="Structure", font=("Helvetica", 11, "bold"))
        self.title.pack(anchor="w")
        self.view = tk.Canvas(right, bg="#eef1f6", highlightthickness=0)
        self.view.pack(fill="both", expand=True, pady=4)
        self.view.bind("<ButtonPress-1>", self.view_press)
        self.view.bind("<B1-Motion>", self.view_drag)
        self.view.bind("<ButtonRelease-1>", self.view_release)
        self.view.bind("<MouseWheel>", lambda e: self.zoom_by(1.1 if e.delta > 0 else 1 / 1.1))
        self.view.bind("<Button-4>", lambda e: self.zoom_by(1.1))       # Linux wheel up
        self.view.bind("<Button-5>", lambda e: self.zoom_by(1 / 1.1))   # Linux wheel down
        self.view.bind("<Configure>", lambda e: self.schedule_view())

        row = ttk.Frame(right)
        row.pack(fill="x")
        ttk.Label(row, text="Model:").pack(side="left")
        for key, label in MODELS:
            ttk.Radiobutton(row, text=label, value=key, variable=self.model,
                            command=self.schedule_view).pack(side="left", padx=(4, 0))
        row2 = ttk.Frame(right)
        row2.pack(fill="x", pady=(4, 0))
        ttk.Checkbutton(row2, text="Atom numbers", variable=self.show_labels,
                        command=self.schedule_view).pack(side="left")
        ttk.Checkbutton(row2, text="Count from 0", variable=self.zero_based,
                        command=self.numbering_changed).pack(side="left", padx=(8, 0))
        ttk.Button(row2, text="Reset view", command=self.reset_view).pack(side="left", padx=(12, 4))
        ttk.Button(row2, text="Clear measurement", command=self.clear_measure).pack(side="left")
        self.measure_lbl = ttk.Label(right, text="", font=("Helvetica", 11, "bold"))
        self.measure_lbl.pack(anchor="w", pady=(4, 0))
        ttk.Label(right, text="Drag to rotate, mouse wheel to zoom. Click 2, 3 or 4 atoms to measure "
                              "a bond, angle or dihedral; click empty space to clear. "
                              "Keys 1–4 change the model, N turns atom numbers on/off.",
                  foreground="#666", wraplength=560).pack(anchor="w")

        keys = (("<Left>", lambda: self.step(-1)), ("<Up>", lambda: self.step(-1)),
                ("<Right>", lambda: self.step(1)), ("<Down>", lambda: self.step(1)),
                ("<Delete>", self.remove_action), ("<BackSpace>", self.remove_action),
                ("<Control-z>", self.undo), ("<Control-o>", self.open_dialog),
                ("<Control-s>", self.save_xyz), ("<Control-q>", self.root.destroy),
                ("<Control-S>", self.save_current), ("<Control-Shift-s>", self.save_current),
                ("n", self.toggle_labels), ("N", self.toggle_labels))
        for seq, fn in keys:
            self.root.bind(seq, lambda e, fn=fn: None if self.typing() else fn())
        for k, (key, _label) in enumerate(MODELS):
            self.root.bind(str(k + 1), lambda e, key=key: None if self.typing() else self.set_model(key))
        # the tables handle arrow keys themselves; keep ours
        for tree in (self.tree, self.vtree):
            for seq, d in (("<Left>", -1), ("<Up>", -1), ("<Right>", 1), ("<Down>", 1)):
                tree.bind(seq, lambda e, d=d: (self.step(d), "break")[1])

    def build_vib_tab(self, vf):
        ttk.Label(vf, text="Click a mode to animate it. Imaginary modes are red. "
                           "Click a peak in the spectrum to pick its mode.",
                  foreground="#666", wraplength=560).pack(anchor="w")
        tf = ttk.Frame(vf)
        tf.pack(fill="x", pady=(4, 0))
        cols = ("no", "freq", "ir")
        self.vtree = ttk.Treeview(tf, columns=cols, show="headings", selectmode="browse", height=8)
        for c, t, w in zip(cols, ("Mode", "Frequency (cm⁻¹)", "IR intensity (km/mol)"), (70, 150, 170)):
            self.vtree.heading(c, text=t)
            self.vtree.column(c, width=w, anchor="e")
        self.vtree.tag_configure("imag", foreground="#c62828")
        sb = ttk.Scrollbar(tf, orient="vertical", command=self.vtree.yview)
        self.vtree.configure(yscrollcommand=sb.set)
        self.vtree.pack(side="left", fill="x", expand=True)
        sb.pack(side="right", fill="y")
        self.vtree.bind("<<TreeviewSelect>>", self.vtree_select)

        self.spec = tk.Canvas(vf, height=230, bg="white", highlightthickness=0)
        self.spec.pack(fill="x", pady=4)
        self.spec.bind("<Button-1>", self.spec_click)
        self.spec.bind("<Configure>", lambda e: self.draw_spec())

        r1 = ttk.Frame(vf)
        r1.pack(fill="x", pady=2)
        self.b_anim = ttk.Button(r1, text="Animate", width=9, command=self.toggle_anim)
        self.b_anim.pack(side="left")
        ttk.Label(r1, text="Amplitude").pack(side="left", padx=(10, 2))
        ttk.Scale(r1, from_=0.1, to=1.0, variable=self.amp, length=110).pack(side="left")
        ttk.Checkbutton(r1, text="Arrows", variable=self.show_arrows,
                        command=self.schedule_view).pack(side="left", padx=(10, 0))
        ttk.Label(r1, text="Peak width (cm⁻¹)").pack(side="left", padx=(12, 2))
        sp = ttk.Spinbox(r1, from_=2, to=200, increment=2, width=5, textvariable=self.fwhm,
                         command=self.draw_spec)
        sp.pack(side="left")
        sp.bind("<Return>", lambda e: self.draw_spec())
        sp.bind("<FocusOut>", lambda e: self.draw_spec())

        box = ttk.LabelFrame(vf, text="Displace along the selected mode (start structures for IRC / follow-up opt)",
                             padding=6)
        box.pack(fill="x", pady=(6, 0))
        ttk.Label(box, text="Step: the atom that moves most moves").pack(side="left")
        ttk.Spinbox(box, from_=0.01, to=1.0, increment=0.01, width=6, format="%.2f",
                    textvariable=self.disp_step).pack(side="left", padx=4)
        ttk.Label(box, text="Å").pack(side="left")
        self.b_disp = ttk.Button(box, text="Save forward + backward .xyz…", command=self.save_displaced)
        self.b_disp.pack(side="left", padx=(12, 0))

    def typing(self):
        w = self.root.focus_get()
        return isinstance(w, (tk.Entry, ttk.Entry, tk.Spinbox))

    def set_model(self, key):
        self.model.set(key)
        self.schedule_view()

    def toggle_labels(self):
        self.show_labels.set(not self.show_labels.get())
        self.schedule_view()

    def numbering_changed(self):
        self.update_scan_heading()
        self.schedule_view()

    def atom_no(self, i):
        return i + (0 if self.zero_based.get() else 1)

    # ---------- loading
    def open_dialog(self):
        path = filedialog.askopenfilename(
            title="Open file", initialdir=os.path.dirname(self.path) or self.here,
            filetypes=[("Gaussian / ORCA / xyz", "*.log *.out *.xyz"), ("All files", "*")])
        if path:
            self.load_path(path)

    def load_path(self, path):
        try:
            with open(path, "r", errors="replace") as fh:
                text = fh.read()
            res = parse_any(os.path.basename(path), text)
        except Exception as exc:  # show any problem to the user, keep the old file open
            messagebox.showerror("Could not open file", "%s\n\n%s" % (os.path.basename(path), exc))
            return False
        self.path, self.name = path, os.path.basename(path)
        self.kind, self.frames, self.scan_coord, self.program = res.kind, res.frames, res.scan_coord, res.program
        self.zero_based.set(res.program == "orca")
        self.removed, self.undo_stack, self.selected, self.current = set(), [], set(), 0
        self.bond_cache = {}
        self.measure = list(self.scan_coord[1]) if self.scan_coord else []
        self.vib, self.vib_sel = res.vib, None
        ext = 1.0
        for atoms in [f.atoms for f in self.frames] + ([self.vib.atoms] if self.vib else []):
            c = centroid(atoms)
            for a in atoms:
                ext = max(ext, length((a[1] - c[0], a[2] - c[1], a[3] - c[2])) + 0.6)
        self.extent = ext
        self.fill_vib()
        self.rot, self.zoom = self.default_rot(), 1.0
        kinds = {"scan": "relaxed scan, %d converged points", "opt": "optimisation, %d steps",
                 "steps": "%d structures", "xyz": "multi-frame .xyz, %d frames", "single": "single structure"}
        k = kinds[self.kind]
        prog = {"gaussian": "Gaussian", "orca": "ORCA", "xyz": "xyz"}.get(self.program, "")
        self.info.config(text="%s — %s %s, %d atoms" % (self.name, prog, k % len(self.frames) if "%d" in k else k,
                                                        len(self.frames[0].atoms)))
        if self.vib:
            n_im = sum(1 for m in self.vib.modes if m.freq < 0)
            self.info.config(text=self.info.cget("text") + " · %d vibrations (%d imaginary)"
                             % (len(self.vib.modes), n_im))
        self.root.title("%s — %s" % (APP, self.name))
        self.fill_tree()
        self.refresh()
        # frequency job or a TS: open the Vibrations tab and animate the first imaginary mode
        if self.vib and (self.kind == "single" or any(m.freq < 0 for m in self.vib.modes)):
            self.nb.select(self.vib_tab)
            self.anim_on = False
            self.toggle_anim()
        else:
            self.nb.select(self.energy_tab)
            self.stop_anim()
        return True

    # ---------- vibrations
    def fill_vib(self):
        self.vtree.delete(*self.vtree.get_children())
        if not self.vib:
            self.nb.tab(self.vib_tab, state="disabled")
            self.draw_spec()
            return
        self.nb.tab(self.vib_tab, state="normal")
        for k, m in enumerate(self.vib.modes):
            self.vtree.insert("", "end", iid=str(k), tags=("imag",) if m.freq < 0 else (),
                              values=(m.label, ("%.2f" % m.freq) + (" (imag.)" if m.freq < 0 else ""),
                                      "%.2f" % m.ir))
        first = next((k for k, m in enumerate(self.vib.modes) if m.freq < 0), 0)
        self.select_mode(first)

    def vib_active(self):
        return self.vib is not None and self.nb.select() == str(self.vib_tab)

    def on_tab(self, _e=None):
        if not hasattr(self, "vtree"):
            return
        if self.vib_active():
            if self.anim_on and self.anim_job is None:
                self.anim_tick()
        self.measure = list(self.scan_coord[1]) if (self.scan_coord and not self.vib_active()) else []
        self.schedule_view()

    def select_mode(self, k):
        if not self.vib or not (0 <= k < len(self.vib.modes)):
            return
        self.vib_sel = k
        self.vib_busy = True
        self.vtree.selection_set(str(k))
        self.vtree.see(str(k))
        self.root.after_idle(self._vib_free)
        self.draw_spec()
        self.schedule_view()

    def _vib_free(self):
        self.vib_busy = False

    def vtree_select(self, _e):
        if self.vib_busy:
            return
        sel = self.vtree.selection()
        if sel:
            self.select_mode(int(sel[0]))
            if not self.anim_on:
                self.toggle_anim()

    def toggle_anim(self):
        if self.anim_on:
            self.stop_anim()
            return
        self.anim_on = True
        self.b_anim.config(text="Stop")
        if self.anim_job is None:
            self.anim_tick()

    def stop_anim(self):
        self.anim_on = False
        self.b_anim.config(text="Animate")
        if self.anim_job is not None:
            self.root.after_cancel(self.anim_job)
            self.anim_job = None
        self.phase = 0.0
        self.schedule_view()

    def anim_tick(self):
        self.anim_job = None
        if not self.anim_on or not self.vib_active():
            return
        self.phase = (self.phase + 2 * math.pi / 30) % (2 * math.pi)
        self.draw_view()
        self.anim_job = self.root.after(40, self.anim_tick)

    def mode(self):
        return self.vib.modes[self.vib_sel] if (self.vib and self.vib_sel is not None) else None

    def number_from(self, var, lo, hi, what):
        try:
            v = float(var.get().replace(",", "."))
        except ValueError:
            v = None
        if v is None or not lo <= v <= hi:
            messagebox.showinfo(APP, "%s must be a number between %g and %g." % (what, lo, hi))
            return None
        return v

    def save_displaced(self, folder=None):
        m = self.mode()
        if m is None:
            return None
        step = self.number_from(self.disp_step, 0.001, 2.0, "The step")
        if step is None:
            return None
        ask = folder is None
        if ask:
            folder = filedialog.askdirectory(title="Folder for the two .xyz files",
                                             initialdir=os.path.dirname(self.path) or ".")
        if not folder:
            return None
        base = os.path.splitext(self.name)[0]
        written = []
        for sign, word in ((1, "fwd"), (-1, "bwd")):
            atoms = displaced(self.vib.atoms, m.vec, sign * step)
            path = os.path.join(folder, "%s_mode%s_%s.xyz" % (base, m.label, word))
            with open(path, "w", newline="\n") as fh:
                fh.write("%d\n%s %s%.3f A along mode %s (%.2f cm-1), largest atom move %.3f A\n"
                         % (len(atoms), self.name, "+" if sign > 0 else "-", step, m.label, m.freq, step))
                for el, x, y, z in atoms:
                    fh.write("%-2s %14.8f %14.8f %14.8f\n" % (el, x, y, z))
            written.append(path)
        self.info.config(text="Saved %s and %s" % tuple(os.path.basename(p) for p in written))
        if ask:
            messagebox.showinfo(APP, "Saved two structures:\n\n%s\n%s\n\nin %s"
                                % (os.path.basename(written[0]), os.path.basename(written[1]), folder))
        return written

    def draw_spec(self):
        c = self.spec
        c.delete("all")
        self.spec_pts = []
        W, H = max(c.winfo_width(), 200), max(c.winfo_height(), 120)
        if not self.vib:
            c.create_text(W / 2, H / 2, text="No frequencies in this file.", fill="#667085")
            return
        real = [m for m in self.vib.modes if m.freq > 0]
        if not real:
            c.create_text(W / 2, H / 2, text="Only imaginary modes — no spectrum.", fill="#667085")
            return
        try:
            fwhm = max(1.0, min(500.0, float(self.fwhm.get().replace(",", "."))))
        except ValueError:
            fwhm = 20.0
        top = max(m.freq for m in real)
        hi = 4000.0 if top <= 3800 else math.ceil((top + 200) / 500.0) * 500
        ml, mr, mt, mb = 50, 14, 14, 36

        def X(nu):                            # high wavenumbers on the left, as usual for IR
            return ml + (hi - nu) / hi * (W - ml - mr)

        curve = ir_curve(real, 0.0, hi, fwhm, max(200, int(W - ml - mr)))
        ymax = max([v for _, v in curve] + [m.ir for m in real]) or 1.0

        def Y(v):
            return H - mb - v / ymax * (H - mt - mb)

        for t in range(0, int(hi) + 1, 500):
            c.create_line(X(t), mt, X(t), H - mb, fill="#eef1f6")
            c.create_text(X(t), H - mb + 12, text=str(t), fill="#667085", font=("Helvetica", 9))
        c.create_line(ml, H - mb, W - mr, H - mb, fill="#9aa3b2")
        c.create_line(ml, mt, ml, H - mb, fill="#9aa3b2")
        c.create_text((ml + W - mr) / 2, H - 8, text="Wavenumber (cm⁻¹)", font=("Helvetica", 10))
        try:
            c.create_text(14, (mt + H - mb) / 2, text="IR intensity", angle=90, font=("Helvetica", 10))
        except tk.TclError:
            pass
        sel = self.mode()
        for m in real:
            x = X(m.freq)
            on = m is sel
            c.create_line(x, Y(0), x, Y(m.ir), fill="#d43c3c" if on else "#9aa3b2", width=3 if on else 1)
            self.spec_pts.append((x, self.vib.modes.index(m)))
        pts = []
        for nu, v in curve:
            pts += [X(nu), Y(v)]
        c.create_line(*pts, fill="#2f6fdf", width=2)
        n_im = len(self.vib.modes) - len(real)
        if n_im:
            c.create_text(W - mr - 4, mt + 2, anchor="ne", fill="#c62828", font=("Helvetica", 9),
                          text="%d imaginary mode%s not shown in the spectrum" % (n_im, "" if n_im == 1 else "s"))

    def spec_click(self, e):
        best = None
        for x, k in self.spec_pts:
            d = abs(x - e.x)
            if d <= 8 and (best is None or d < best[0]):
                best = (d, k)
        if best:
            self.select_mode(best[1])
            if not self.anim_on:
                self.toggle_anim()

    def x_label(self):
        return self.LABELS.get(self.kind, "Point")

    def visible(self):
        return [i for i in range(len(self.frames)) if i not in self.removed]

    def energy_ref(self):
        es = [self.frames[i].energy for i in self.visible() if self.frames[i].energy is not None]
        return min(es) if es else None

    # ---------- removing points
    def remove_points(self, idx):
        idx = [i for i in idx if i not in self.removed]
        if not idx:
            return
        if len(idx) >= len(self.visible()):
            messagebox.showinfo(APP, "At least one point has to stay.")
            return
        self.removed.update(idx)
        self.undo_stack.append(idx)
        self.selected.clear()
        if self.current in self.removed:
            vis = self.visible()
            self.current = min(vis, key=lambda i: abs(i - self.current))
        self.refresh()

    def remove_action(self):
        if len(self.frames) > 1 and not self.vib_active():
            self.remove_points(sorted(self.selected) if self.selected else [self.current])

    def undo(self):
        if self.undo_stack:
            for i in self.undo_stack.pop():
                self.removed.discard(i)
            self.refresh()

    def restore_all(self):
        self.removed.clear()
        self.undo_stack = []
        self.refresh()

    def bring_back(self, i):
        self.removed.discard(i)
        self.undo_stack = [[x for x in l if x != i] for l in self.undo_stack]
        self.undo_stack = [l for l in self.undo_stack if l]
        self.refresh()

    def step(self, d):
        if self.vib_active():
            if self.vib_sel is not None:
                self.select_mode(max(0, min(len(self.vib.modes) - 1, self.vib_sel + d)))
            return
        vis = self.visible()
        if not vis:
            return
        k = vis.index(self.current) if self.current in vis else 0
        self.current = vis[max(0, min(len(vis) - 1, k + d))]
        self.selected.clear()
        self.refresh()

    # ---------- saving
    def comment(self, i):
        f = self.frames[i]
        e = "Energy: %.10f Eh | " % f.energy if f.energy is not None else ""
        return "%s%s %d | %s" % (e, self.x_label().lower(), i + 1, self.name)

    def xyz_lines(self):
        out = []
        for i in self.visible():
            f = self.frames[i]
            out.append(str(len(f.atoms)))
            out.append(self.comment(i))
            for el, x, y, z in f.atoms:
                out.append("%-2s %14.8f %14.8f %14.8f" % (el, x, y, z))
        return out

    def save_xyz(self, path=None):
        if not self.frames:
            return
        if path is None:
            base = os.path.splitext(self.name)[0]
            path = filedialog.asksaveasfilename(
                title="Save", initialdir=os.path.dirname(self.path) or ".",
                initialfile=base + "-cleaned.xyz", defaultextension=".xyz", filetypes=[(".xyz", "*.xyz")])
        if path:
            with open(path, "w", newline="\n") as fh:
                fh.write("\n".join(self.xyz_lines()) + "\n")
            self.info.config(text="Saved %d points to %s" % (len(self.visible()), os.path.basename(path)))

    def save_current(self, path=None):
        """Save only the structure shown on the right (resting structure in the Vibrations tab)."""
        if not self.frames:
            return None
        base = os.path.splitext(self.name)[0]
        if self.vib_active():
            atoms, tag = self.vib.atoms, "freq-geom"
            comment = "%s | structure of the frequency calculation" % self.name
        else:
            f = self.frames[self.current]
            atoms, tag = f.atoms, "%s%d" % ("point" if self.kind == "scan" else "step", self.current + 1)
            comment = self.comment(self.current)
            if len(self.frames) == 1:
                tag = "geom"
        if path is None:
            path = filedialog.asksaveasfilename(
                title="Save this structure", initialdir=os.path.dirname(self.path) or ".",
                initialfile="%s_%s.xyz" % (base, tag), defaultextension=".xyz", filetypes=[(".xyz", "*.xyz")])
        if not path:
            return None
        with open(path, "w", newline="\n") as fh:
            fh.write("%d\n%s\n" % (len(atoms), comment))
            for el, x, y, z in atoms:
                fh.write("%-2s %14.8f %14.8f %14.8f\n" % (el, x, y, z))
        self.info.config(text="Saved one structure to %s" % os.path.basename(path))
        return path

    # ---------- table
    def update_scan_heading(self):
        if self.scan_coord:
            t, at = self.scan_coord
            name = {"B": "Bond", "A": "Angle", "D": "Dihedral", "L": "Linear"}.get(t, t)
            self.tree.heading("coord", text="%s %s" % (name, "-".join(str(self.atom_no(a)) for a in at)))

    def fill_tree(self):
        self.tree.delete(*self.tree.get_children())
        self.tree["displaycolumns"] = ("pt", "coord", "e", "de", "status") if self.scan_coord \
            else ("pt", "e", "de", "status")
        self.update_scan_heading()
        for i in range(len(self.frames)):
            self.tree.insert("", "end", iid=str(i), values=("",) * 5)

    def update_tree(self):
        ref = self.energy_ref()
        for i, f in enumerate(self.frames):
            coord = ""
            if self.scan_coord:
                mv = measure_value(f.atoms, self.scan_coord[1])
                if mv:
                    coord = ("%.3f " if mv[2] == "Å" else "%.1f ") % mv[1] + mv[2]
            rem = i in self.removed
            if rem:
                de = "removed"
            elif f.energy is not None and ref is not None:
                de = "%.2f" % ((f.energy - ref) * HARTREE_TO_KCAL)
            else:
                de = "—"
            e = "%.8f" % f.energy if f.energy is not None else "—"
            self.tree.item(str(i), values=(i + 1, coord, e, de, f.status), tags=("removed",) if rem else ())
        self.tree_busy = True
        sel = [str(i) for i in sorted(self.selected)] or [str(self.current)]
        self.tree.selection_set(sel)
        self.tree.see(str(self.current))
        self.root.after_idle(self._tree_free)

    def _tree_free(self):
        self.tree_busy = False

    def tree_select(self, _e):
        if self.tree_busy:
            return
        sel = [int(s) for s in self.tree.selection() if int(s) not in self.removed]
        if not sel:
            return
        self.selected = set(sel) if len(sel) > 1 else set()
        self.current = sel[-1]
        self.refresh(tree=False)

    def tree_double(self, e):
        row = self.tree.identify_row(e.y)
        if row and int(row) in self.removed:
            self.bring_back(int(row))

    # ---------- plot
    def draw_plot(self):
        c = self.plot
        c.delete("all")
        W, H = max(c.winfo_width(), 200), max(c.winfo_height(), 150)
        self.plot_pts = []
        ref = self.energy_ref()
        vis = [i for i in self.visible() if self.frames[i].energy is not None]
        if len(self.frames) < 2 or not vis:
            msg = "Open a file with Ctrl+O" if not self.frames else (
                "Only one structure — nothing to plot." if len(self.frames) < 2 else "No energies to plot.")
            c.create_text(W / 2, H / 2, text=msg, fill="#667085")
            return
        ml, mr, mt, mb = 62, 16, 14, 40
        rel = {i: (self.frames[i].energy - ref) * HARTREE_TO_KCAL for i in vis}
        ymax = max(rel.values()) or 1.0
        pad = ymax * 0.06
        n = len(self.frames)

        def X(v):
            return ml + (v - 1) / float(max(1, n - 1)) * (W - ml - mr)

        def Y(v):
            return H - mb - (v + pad) / (ymax + 2 * pad) * (H - mt - mb)

        for t in nice_ticks(0, ymax, 5):
            c.create_line(ml, Y(t), W - mr, Y(t), fill="#eef1f6")
            c.create_text(ml - 6, Y(t), text=("%g" % round(t, 3)), anchor="e", fill="#667085", font=("Helvetica", 9))
        for t in nice_ticks(1, n, min(10, max(1, n - 1))):
            if float(t).is_integer():
                c.create_text(X(t), H - mb + 12, text=str(int(t)), fill="#667085", font=("Helvetica", 9))
        c.create_line(ml, H - mb, W - mr, H - mb, fill="#9aa3b2")
        c.create_line(ml, mt, ml, H - mb, fill="#9aa3b2")
        c.create_text((ml + W - mr) / 2, H - 8, text=self.x_label(), font=("Helvetica", 10))
        try:
            c.create_text(14, (mt + H - mb) / 2, text="ΔE (kcal/mol)", angle=90, font=("Helvetica", 10))
        except tk.TclError:
            c.create_text(30, mt + 4, text="ΔE (kcal/mol)", anchor="w", font=("Helvetica", 10))
        if len(vis) > 1:
            pts = []
            for i in vis:
                pts += [X(i + 1), Y(rel[i])]
            c.create_line(*pts, fill="#2f6fdf", width=2)
        for i in vis:
            x, y = X(i + 1), Y(rel[i])
            self.plot_pts.append((x, y, i))
            cur, sel = i == self.current, i in self.selected
            bad = self.frames[i].status == "not converged"
            if cur:
                c.create_oval(x - 9, y - 9, x + 9, y + 9, outline="#f5c400", width=3)
            r = 6 if cur else 4.5
            c.create_oval(x - r, y - r, x + r, y + r,
                          fill="#d43c3c" if sel else ("white" if bad else "#2f6fdf"),
                          outline="#e08a00" if bad else "white", width=2 if bad else 1)

    def plot_click(self, e, ctrl):
        best = None
        for x, y, i in self.plot_pts:
            d = math.hypot(x - e.x, y - e.y)
            if d <= 12 and (best is None or d < best[0]):
                best = (d, i)
        if not best:
            return
        i = best[1]
        if ctrl:
            self.selected.symmetric_difference_update({i})
        else:
            self.selected.clear()
        self.current = i
        self.refresh()

    # ---------- 3D view
    @staticmethod
    def default_rot():
        return mat_mul(rot_x(-0.5), rot_y(0.6))

    def reset_view(self):
        self.rot, self.zoom = self.default_rot(), 1.0
        self.schedule_view()

    def zoom_by(self, f):
        self.zoom = max(0.2, min(8.0, self.zoom * f))
        self.schedule_view()

    def schedule_view(self):
        if not self.redraw_pending:
            self.redraw_pending = True
            self.root.after_idle(self.draw_view)

    def draw_view(self):
        self.redraw_pending = False
        c = self.view
        c.delete("all")
        W, H = max(c.winfo_width(), 100), max(c.winfo_height(), 100)
        if not self.frames:
            self.projected = []
            return
        model = self.model.get()
        mode = self.mode() if self.vib_active() else None
        if self.vib_active():
            eq = self.vib.atoms
            cx, cy, cz = centroid(eq)
            atoms = eq
            if mode is not None and self.anim_on:
                atoms = displaced(eq, mode.vec, self.amp.get() * math.sin(self.phase))
            cache_key = "vib"
            self.shown_atoms = eq             # measure on the resting structure
            if mode is not None:
                self.title.config(text="Vibration — mode %s: %.2f cm⁻¹%s" % (
                    mode.label, mode.freq, " (imaginary)" if mode.freq < 0 else ""))
        else:
            atoms = self.frames[self.current].atoms
            cx, cy, cz = centroid(atoms)
            cache_key = self.current
            self.shown_atoms = atoms
            f = self.frames[self.current]
            self.title.config(text="Structure — %s %d%s" % (self.x_label().lower(), self.current + 1,
                                                           " (%s)" % f.status if f.status else ""))
        R = self.rot
        extent = self.extent + (1.3 if model == "space" else 0.0)
        s = min(W, H) / (2 * extent) * self.zoom
        proj = []
        for i, (el, x, y, z) in enumerate(atoms):
            v = (x - cx, y - cy, z - cz)
            p = [R[k][0] * v[0] + R[k][1] * v[1] + R[k][2] * v[2] for k in range(3)]
            proj.append({"i": i, "el": el, "sx": W / 2 + p[0] * s, "sy": H / 2 - p[1] * s, "z": p[2],
                         "r": atom_radius(el, model) * s})
        self.projected = proj
        bonds = self.bond_cache.get(cache_key)
        if bonds is None:                     # bonds from the resting structure, so they don't flicker
            bonds = self.bond_cache[cache_key] = bonds_of(self.vib.atoms if cache_key == "vib" else atoms)

        items = []
        if model != "space":
            for i, j in bonds:
                A, B = proj[i], proj[j]
                dx, dy = B["sx"] - A["sx"], B["sy"] - A["sy"]
                L = math.hypot(dx, dy)
                if model == "ball":
                    if L <= A["r"] + B["r"] + 1:
                        continue              # bond points at the viewer, hidden by the atoms
                    ux, uy = dx / L, dy / L
                    ax, ay = A["sx"] + ux * (A["r"] - 0.5), A["sy"] + uy * (A["r"] - 0.5)
                    bx, by = B["sx"] - ux * (B["r"] - 0.5), B["sy"] - uy * (B["r"] - 0.5)
                else:                         # sticks and wireframe start at the atom centre
                    if L < 1:
                        continue
                    ax, ay, bx, by = A["sx"], A["sy"], B["sx"], B["sy"]
                mx, my = (ax + bx) / 2, (ay + by) / 2
                items.append(((A["z"] * 3 + B["z"]) / 4 - 0.01, "b", (ax, ay, mx, my, color_of(A["el"]))))
                items.append(((B["z"] * 3 + A["z"]) / 4 - 0.01, "b", (bx, by, mx, my, color_of(B["el"]))))
        if model in ("ball", "space"):
            for p in proj:
                items.append((p["z"], "a", p))
        elif model == "stick":
            for p in proj:                    # small caps on top of their own bonds
                items.append((p["z"] + 0.05, "a", p))
        else:                                 # wireframe: mark atoms without bonds with a cross
            bonded = {k for b in bonds for k in b}
            for p in proj:
                if p["i"] not in bonded:
                    items.append((p["z"], "x", p))
        items.sort(key=lambda t: t[0])

        bw = {"ball": max(2, 0.16 * s), "stick": max(3, 2 * STICK_R * s), "wire": 2}.get(model, 2)
        for _z, kind, d in items:
            if kind == "b":
                x1, y1, x2, y2, col = d
                if model == "wire":
                    c.create_line(x1, y1, x2, y2, fill=shade(col, -0.3), width=bw, capstyle="round")
                else:
                    c.create_line(x1, y1, x2, y2, fill="#555c68", width=bw + 2, capstyle="butt")
                    c.create_line(x1, y1, x2, y2, fill=shade(col, -0.05), width=bw, capstyle="butt")
            elif kind == "x":
                x, y, col = d["sx"], d["sy"], shade(color_of(d["el"]), -0.3)
                c.create_line(x - 5, y - 5, x + 5, y + 5, fill=col, width=2)
                c.create_line(x - 5, y + 5, x + 5, y - 5, fill=col, width=2)
            elif model == "stick":
                x, y, r = d["sx"], d["sy"], bw / 2 + 1
                c.create_oval(x - r, y - r, x + r, y + r, fill=shade(color_of(d["el"]), -0.05), outline="#555c68")
            else:
                x, y, r, col = d["sx"], d["sy"], d["r"], color_of(d["el"])
                c.create_oval(x - r, y - r, x + r, y + r, fill=shade(col, -0.15), outline="#3a4250")
                c.create_oval(x - r * 0.85, y - r * 0.9, x + r * 0.7, y + r * 0.65, fill=col, outline="")
                c.create_oval(x - r * 0.55, y - r * 0.6, x - r * 0.05, y - r * 0.15, fill=shade(col, 0.6), outline="")
        if mode is not None and self.show_arrows.get():   # displacement arrows
            mm = max_move(mode.vec) or 1.0
            for p, v in zip(proj, mode.vec):
                q = [R[k][0] * v[0] + R[k][1] * v[1] + R[k][2] * v[2] for k in range(3)]
                dx, dy = q[0] / mm * 1.0 * s, -q[1] / mm * 1.0 * s
                if math.hypot(dx, dy) < 4:
                    continue
                c.create_line(p["sx"], p["sy"], p["sx"] + dx, p["sy"] + dy, fill="#d4380d", width=2.5,
                              arrow="last", arrowshape=(10, 12, 4))

        def hidden(p):                        # spacefill: centre covered by a closer atom
            return model == "space" and any(
                q["z"] > p["z"] and math.hypot(q["sx"] - p["sx"], q["sy"] - p["sy"]) < q["r"] * 0.8 for q in proj)

        for p in proj:                        # picked atoms
            if p["i"] in self.measure and not hidden(p):
                r = (min(p["r"], 14) if model == "space" else max(p["r"], bw / 2 + 1, 6)) + 3
                c.create_oval(p["sx"] - r, p["sy"] - r, p["sx"] + r, p["sy"] + r, outline="#f5c400", width=3)
        if self.show_labels.get():           # numbers last, so nothing covers them
            for p in sorted(proj, key=lambda q: q["z"]):
                if hidden(p):
                    continue
                size = int(max(9, min(14, p["r"] * 0.8)))
                f = ("Helvetica", size, "bold")
                t = str(self.atom_no(p["i"]))
                for ox, oy in ((-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)):
                    c.create_text(p["sx"] + ox, p["sy"] + oy, text=t, fill="white", font=f)
                c.create_text(p["sx"], p["sy"], text=t, fill="#111111", font=f)
        self.update_measure()

    def view_press(self, e):
        self.drag = [e.x, e.y, 0]

    def view_drag(self, e):
        if not self.drag:
            return
        dx, dy = e.x - self.drag[0], e.y - self.drag[1]
        self.drag = [e.x, e.y, self.drag[2] + abs(dx) + abs(dy)]
        self.rot = mat_mul(mat_mul(rot_y(dx * 0.01), rot_x(dy * 0.01)), self.rot)
        self.schedule_view()

    def view_release(self, e):
        if self.drag and self.drag[2] < 5:
            self.pick_at(e.x, e.y)
        self.drag = None

    def pick_at(self, x, y):
        hit = None
        for p in self.projected:
            if math.hypot(p["sx"] - x, p["sy"] - y) <= max(p["r"], 7) + 2 and (hit is None or p["z"] > hit["z"]):
                hit = p
        if hit is None:
            self.measure = []
        elif hit["i"] in self.measure:
            self.measure.remove(hit["i"])
        else:
            self.measure = [hit["i"]] if len(self.measure) >= 4 else self.measure + [hit["i"]]
        self.schedule_view()

    def clear_measure(self):
        self.measure = []
        self.schedule_view()

    def update_measure(self):
        txt = ""
        atoms = getattr(self, "shown_atoms", None)
        if atoms:
            if len(self.measure) >= 2:
                mv = measure_value(atoms, self.measure)
                if mv:
                    names = "–".join("%s%d" % (atoms[i][0], self.atom_no(i)) for i in self.measure)
                    txt = "%s %s: %s %s" % (mv[0], names, ("%.3f" if mv[2] == "Å" else "%.1f") % mv[1], mv[2])
            elif len(self.measure) == 1:
                i = self.measure[0]
                txt = "Picked %s%d — click one more atom." % (atoms[i][0], self.atom_no(i))
        self.measure_lbl.config(text=txt)

    # ---------- refresh everything
    def refresh(self, tree=True):
        if not self.frames:
            return
        self.b_remove.config(text=("Remove %d marked points" % len(self.selected)) if self.selected
                             else "Remove point %d" % (self.current + 1))
        self.b_undo.state(["!disabled"] if self.undo_stack else ["disabled"])
        self.b_restore.state(["!disabled"] if self.removed else ["disabled"])
        self.draw_plot()
        if tree:
            self.update_tree()
        self.schedule_view()


def main():
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(__doc__)
        return
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        sys.exit("Cannot open a window (%s).\nAre you in a graphical session?" % exc)
    app = App(root)
    if len(sys.argv) > 1:
        root.after(50, lambda: app.load_path(sys.argv[1]))
    else:
        sample = os.path.join(app.here, "sample-scan.log")
        if os.path.exists(sample):
            root.after(50, lambda: app.load_path(sample))
    root.mainloop()


if __name__ == "__main__":
    main()
