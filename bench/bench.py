import shutil, subprocess, sys, time, os, glob

BENCH = os.path.dirname(os.path.abspath(__file__))
os.chdir(BENCH)

def restore(model):
    for f in glob.glob(os.path.join('meshes', model + '.*')):
        shutil.copy(f, os.path.basename(f))

def run(exe, model, reps):
    times = []
    for _ in range(reps):
        restore(model)
        t0 = time.perf_counter()
        r = subprocess.run([exe, model], capture_output=True, text=True)
        t1 = time.perf_counter()
        if r.returncode != 0:
            print(f'  ERROR rc={r.returncode}: {r.stdout[-300:]} {r.stderr[-300:]}')
            return None
        times.append(t1 - t0)
    return times

if __name__ == '__main__':
    exe = sys.argv[1]
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    models = sys.argv[3].split(',') if len(sys.argv) > 3 else ['temp', 'tq', 'age', 'temp_dense', 'tq_dense', 'age_dense']
    print(f'=== {exe} ===')
    for m in models:
        ts = run(exe, m, reps)
        if ts:
            print(f'{m:12s} best={min(ts):7.3f}s  all=' + ' '.join(f'{t:.3f}' for t in ts))
