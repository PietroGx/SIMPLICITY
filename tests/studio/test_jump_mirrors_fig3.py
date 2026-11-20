"""export_transmissions' per-event jumps must be fig3's violin data, event for event.

One simulation per case: the exporter's per-day jump table for that transmitter
type, expanded to one value per event, against
fig3_preprocess_data.get_fig2_divergence_data restricted to the same seed. Reads
only; writes nothing.

    python tests/studio/test_jump_mirrors_fig3.py --exp-num 2
"""
import argparse, base64, os, sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tests"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "nature_plots"))
import export_transmissions as ex
import fig3_preprocess_data as f3
import simplicity.dir_manager as dm
from _scenarios import EXP_NAME, one_sod


def decode(rec, t, nd):
    j = rec["jumps"][t]
    H = np.frombuffer(base64.b64decode(j["h"]), np.uint16).reshape(nd, j["D"]) if j["D"] else np.zeros((nd, 0))
    return np.repeat(np.arange(H.shape[1]), H.sum(0).astype(int))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-num", type=int, required=True)
    a = ap.parse_args()
    # fig3 reads long shedders' jumps everywhere but control, standard hosts' in control
    cases = [("SOT", "long_shedder", ("ll", "ls")), ("HIV_high", "long_shedder", ("ll", "ls")),
             ("control", "standard", ("sl", "ss"))]
    bad = 0
    for sc, t, cells in cases:
        ssod = dm.get_seeded_simulation_output_dirs(one_sod(EXP_NAME, sc, a.exp_num))[0]
        seed = int(dm.get_seed_from_SSOD(str(ssod)))
        _, _, rec, err = ex._one((sc, ssod))
        assert rec, err
        ours = decode(rec, t, rec["n_days"])
        theirs = sorted(int(v) for v in f3.get_fig2_divergence_data(a.exp_num, sc, {seed}))
        events = int(sum(np.frombuffer(base64.b64decode(rec[c]), np.uint16).sum() for c in cells))
        same = ours.tolist() == theirs and len(ours) == events
        bad += not same
        print(f"{'ok  ' if same else 'FAIL'} {sc} {t} seed {seed}: {len(ours)} jumps "
              f"(fig3 {len(theirs)}, events {events}), median {np.median(ours) if len(ours) else '-'}, "
              f"max {ours.max() if len(ours) else '-'}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
