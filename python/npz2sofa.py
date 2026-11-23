#!/usr/bin/env python3
# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
"""
npz2sofa.py
Converts the NPZ file written by hrtf2hrir.py into a SOFA file.

Convention: SimpleFreeFieldHRIR (AES69-2022)
  Data.IR shape: (M, R, N)  M = directions, R = 2 ears, N = samples
  R = 0: left ear,  R = 1: right ear

Requires: pip install sofar

Usage:
  python npz2sofa.py hrir.npz output.sofa [options]

Metadata (title, author, organisation, licence, database) come from
sofa_meta.txt in the current directory, else python/sofa_meta.txt; the
options override them.

Spatial convention (from hrtf2hrir.py):
  az = atan2(y, x), x front, z up, counter-clockwise
  -> SOFA: az in [0, 360[, el in [-90, 90], r in metres
"""

import argparse
import numpy as np
from pathlib import Path

try:
    import sofar as sf
except ImportError:
    raise ImportError("sofar is not installed: pip install sofar")


META_FILE = 'sofa_meta.txt'
META_DEFAULTS = {
    'title':   'HRTF MLFMM',
    'author':  '',
    'org':     'CNAM',
    'license': 'CC BY 4.0',
    'db':      '',
    'subject': '',
    'comment': '',
}


def load_meta_file(path=None):
    """
    Lit sofa_meta.txt : lignes 'clé = valeur'. Sans --meta : celui du répertoire
    courant s'il existe, sinon celui du répertoire du script (python/sofa_meta.txt).
    Lignes vides et commentaires (# ...) ignorés.
    Retourne un dict avec les valeurs trouvées.
    """
    if path:
        meta_path = Path(path)
    else:
        meta_path = Path(META_FILE)
        if not meta_path.exists():
            meta_path = Path(__file__).resolve().parent / META_FILE
    if not meta_path.exists():
        print(f"No {META_FILE} found: default metadata.")
        return {}
    print(f"Metadata read from: {meta_path}")
    meta = {}
    for line in meta_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '=' in line:
            key, _, val = line.partition('=')
            meta[key.strip().lower()] = val.strip()
    return meta


def load_npz(npz_file):
    """Charge le NPZ produit par hrtf2hrir.py."""
    data = np.load(npz_file)
    required = {'hrir_L', 'hrir_R', 'az', 'el', 'r', 'fe'}
    missing = required - set(data.files)
    if missing:
        raise ValueError(f"{npz_file}: missing fields: {missing}")
    return data


def build_sofa(hrir_L, hrir_R, az, el, r, fe, equalized, meta):
    """
    Construit un objet SOFA SimpleFreeFieldHRIR.

    hrir_L, hrir_R : (N_out, N_ext)
    az, el         : azimut et élévation en degrés
    r              : distance en mètres
    fe             : fréquence d'échantillonnage
    equalized      : booléen, égalisation champ diffus appliquée
    meta           : dict de métadonnées GLOBAL_*
    """
    N_out, N_ext = hrir_L.shape

    # Data.IR : (M, R, N) — M=directions, R=2, N=samples
    # R=0 gauche, R=1 droite
    data_ir = np.stack([hrir_L.T, hrir_R.T], axis=1).astype(np.float64)

    # SourcePosition SOFA : azimut [0, 360[, élévation [-90, 90], rayon (m)
    src_pos = np.column_stack([az % 360, el, r])  # (N_ext, 3)

    sofa = sf.Sofa('SimpleFreeFieldHRIR')

    # Données
    sofa.Data_IR           = data_ir
    sofa.Data_SamplingRate = float(fe)

    # Positions source
    sofa.SourcePosition       = src_pos
    sofa.SourcePosition_Type  = 'spherical'
    sofa.SourcePosition_Units = 'degree, degree, metre'

    # Position/orientation auditeur (origine, regardant +x, haut = +z)
    sofa.ListenerPosition = np.array([[0.0, 0.0, 0.0]])
    sofa.ListenerView     = np.array([[1.0, 0.0, 0.0]])
    sofa.ListenerUp       = np.array([[0.0, 0.0, 1.0]])

    # Métadonnées globales standard
    sofa.GLOBAL_Title         = meta.get('title',   'HRTF MLFMM')
    sofa.GLOBAL_AuthorContact = meta.get('author',  '')
    sofa.GLOBAL_Organization  = meta.get('org',     '')
    sofa.GLOBAL_License       = meta.get('license', 'CC BY 4.0')
    sofa.GLOBAL_Comment       = meta.get('comment', '')
    sofa.GLOBAL_DatabaseName      = meta.get('db',      '')
    sofa.GLOBAL_ListenerShortName = meta.get('subject', '')

    # Attributs personnalisés : traçabilité du traitement MLFMM
    sofa.add_attribute('GLOBAL_DiffuseFieldEqualized', 'yes' if equalized else 'no')
    sofa.GLOBAL_ApplicationName    = 'MLFMM2HRTF'
    sofa.GLOBAL_ApplicationVersion = '1.0'

    return sofa


def main():
    parser = argparse.ArgumentParser(
        description='NPZ (hrtf2hrir.py) -> SOFA file, SimpleFreeFieldHRIR.')
    parser.add_argument('npz_file',  help='input NPZ file (hrtf2hrir.py)')
    parser.add_argument('sofa_file', help='output SOFA file (.sofa)')
    parser.add_argument('--meta',    default=None,
                        help=f'metadata file (default: {META_FILE} in the current directory, '
                             f'else python/{META_FILE})')
    parser.add_argument('--title',   default=None,
                        help='title (GLOBAL_Title), overrides sofa_meta.txt')
    parser.add_argument('--author',  default=None,
                        help='author contact, overrides sofa_meta.txt')
    parser.add_argument('--org',     default=None,
                        help='organisation, overrides sofa_meta.txt')
    parser.add_argument('--license', default=None,
                        help='licence, overrides sofa_meta.txt')
    parser.add_argument('--comment', default=None,
                        help='free comment, overrides sofa_meta.txt')
    parser.add_argument('--db',      default=None,
                        help='database name, overrides sofa_meta.txt')
    parser.add_argument('--subject', default=None,
                        help='subject id, e.g. "S01", overrides sofa_meta.txt')
    args = parser.parse_args()

    # --- Métadonnées : défauts → fichier → arguments CLI ---
    file_meta = load_meta_file(args.meta)

    def resolve(key):
        """Priorité : argument CLI > sofa_meta.txt > défaut."""
        cli_val = getattr(args, key.replace('-', '_'))
        if cli_val is not None:
            return cli_val
        return file_meta.get(key, META_DEFAULTS.get(key, ''))

    # --- Chargement NPZ ---
    data      = load_npz(args.npz_file)
    hrir_L    = data['hrir_L']          # (N_out, N_ext)
    hrir_R    = data['hrir_R']
    az        = data['az']              # (N_ext,) degrés
    el        = data['el']
    r         = data['r']               # (N_ext,) mètres
    fe        = float(data['fe'])
    equalized = bool(data['equalized'])

    N_out, N_ext = hrir_L.shape
    print(f"NPZ read: {args.npz_file}")
    print(f"  N_ext={N_ext} directions | N_out={N_out} samples | fe={fe} Hz")
    print(f"  diffuse-field equalised: {equalized}")

    # --- Construction SOFA ---
    meta = {k: resolve(k) for k in META_DEFAULTS}
    sofa = build_sofa(hrir_L, hrir_R, az, el, r, fe, equalized, meta)

    # --- Écriture ---
    out = Path(args.sofa_file)
    sf.write_sofa(str(out), sofa)
    print(f"SOFA written: {out}")
    print(f"  convention: SimpleFreeFieldHRIR (AES69-2022)")
    print(f"  Data.IR    : {sofa.Data_IR.shape}  (M={N_ext}, R=2, N={N_out})")
    print(f"  licence   : {meta['license']}")
    print(f"  DiffuseFieldEqualized : {'yes' if equalized else 'no'}")


if __name__ == '__main__':
    main()
