"""Generate bounce-fold trunk skeletons with a wall-clock budget.
Run repeatedly; it saves skels_bounce.pkl when it completes.
Usage: python3 gen_skeletons.py [want]
"""
import pickle, sys, time
import witness_tree as wt

want = int(sys.argv[1]) if len(sys.argv) > 1 else 40
t0 = time.time()
sk = wt.build_trunk_skeletons(want=want, filter_live=False)
print(f"skeletons: {len(sk)}  ({time.time()-t0:.0f}s)  "
      f"cache={len(wt._CDP_CACHE)}")
if sk:
    pickle.dump(sk, open("skels_bounce.pkl", "wb"))
    nb = [s.get("nbounce", 0) for s in sk]
    print("bounce counts:", {b: nb.count(b) for b in sorted(set(nb))})
    print("saved skels_bounce.pkl")
