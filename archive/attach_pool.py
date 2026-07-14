"""Parallel attachment over a skeleton pool. Writes progress to
attach_state.json and the first success to witness_tree.json."""
import json, os, pickle, sys, time
from multiprocessing import Pool

import witness_tree as wt

def try_one(i):
    sk = pickle.load(open(POOL, 'rb'))[i]
    t0 = time.time()
    sol = wt.attach_branches(sk)
    return i, sol, time.time() - t0

POOL = sys.argv[1] if len(sys.argv) > 1 else 'skels9.pkl'

if __name__ == '__main__':
    start = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    stop = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    n = len(pickle.load(open(POOL, 'rb')))
    idxs = list(range(start, min(stop, n)))
    with Pool(4) as p:
        for i, sol, dt in p.imap_unordered(try_one, idxs):
            print(f'skeleton {i}:', 'SUCCESS' if sol else 'failed',
                  f'({dt:.0f}s)', flush=True)
            if sol:
                json.dump(sol, open('witness_tree.json', 'w'),
                          indent=1, default=list)
                print('SAVED', flush=True)
                p.terminate()
                break
