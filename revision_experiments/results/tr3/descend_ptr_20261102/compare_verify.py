"""Compare the verification runs of the 2026-11-02 descending search with the published result files.
V2: real_out/real_hd.csv vs results/tr3/real_hd.csv (variant hd, UN-CCDs): n, d, ari, sil, k_hat as stored strings.
V3: sim_out/sim_replicates.csv vs results/tr3/sim_replicates.csv (variant published, UN-CCDs, d 10/20):
    per (dist, d, n, n_clusters, rep): ari, sil, k_hat, data_fp; plus median seconds per cell, old and new.
Writes compare_summary.json (all numbers used later) and prints a report.
Usage: python compare_verify.py <results/tr3> <desc_swap dir>"""
import csv, json, math, os, statistics, sys
from collections import defaultdict
RES, S = sys.argv[1], sys.argv[2]
out = {}


def rows(p):
    return list(csv.DictReader(open(p, newline='', encoding='utf-8')))


# ---- V2 real data
pub = [r for r in rows(os.path.join(RES, 'real_hd.csv')) if r['variant'] == 'hd' and r['method'] == 'UN-CCDs']
new = rows(os.path.join(S, 'real_out', 'real_hd.csv'))
pk = defaultdict(list)
for r in pub:
    pk[r['dataset']].append(r)
v2 = []
for r in new:
    p = pk[r['dataset']]
    assert len(p) == 1, (r['dataset'], len(p))
    p = p[0]
    diffs = [c for c in ('n', 'd', 'ari', 'sil', 'k_hat') if p[c] != r[c]]
    v2.append({'dataset': r['dataset'], 'n': int(r['n']), 'd': int(r['d']), 'equal': not diffs, 'diffs': diffs,
               'sec_old': float(p['seconds']), 'sec_new': float(r['seconds'])})
out['v2_real'] = v2
print('V2 real data (UN-CCDs, single process):')
for x in v2:
    print(f"  {x['dataset']:13s} n={x['n']:5d} d={x['d']:2d} equal={x['equal']} {x['diffs'] or ''} "
          f"seconds {x['sec_old']:.1f} -> {x['sec_new']:.1f}")
print(f"  {sum(x['equal'] for x in v2)} of {len(v2)} equal (6 expected)")

# ---- V3 simulations
if not os.path.exists(os.path.join(S, 'sim_out', 'sim_replicates.csv')):
    print('\nV3: no sim_out/sim_replicates.csv yet'); sys.exit(0)
key =lambda r: (r['dist'], r['d'], r['n'], r['n_clusters'], r['rep'])
pubs = {key(r): r for r in rows(os.path.join(RES, 'sim_replicates.csv'))
        if r['variant'] == 'published' and r['method'] == 'UN-CCDs' and r['d'] in ('10', '20')}
news = rows(os.path.join(S, 'sim_out', 'sim_replicates.csv'))
cnt = defaultdict(lambda: [0, 0])
med = defaultdict(dict)
sec_old = defaultdict(list); sec_new = defaultdict(list)
missing = 0
for r in news:
    p = pubs.get(key(r))
    if p is None:
        missing += 1; continue
    c = (r['dist'], r['d'], r['n'], r['n_clusters'])
    cnt[c][0] += 1
    if all(p[f] == r[f] for f in ('ari', 'sil', 'k_hat', 'data_fp')):
        cnt[c][1] += 1
    sec_old[c].append(float(p['seconds'])); sec_new[c].append(float(r['seconds']))
cells = sorted(cnt, key=lambda c: (int(c[2]), c[0], int(c[1]), int(c[3])))
v3 = []
for c in cells:
    v3.append({'dist': c[0], 'd': int(c[1]), 'n': int(c[2]), 'k': int(c[3]), 'reps': cnt[c][0], 'equal': cnt[c][1],
               'median_sec_old': statistics.median(sec_old[c]), 'median_sec_new': statistics.median(sec_new[c])})
out['v3_sim'] = v3
out['v3_new_rows'] = len(news); out['v3_rows_without_published_match'] = missing
print(f'\nV3 simulations: {len(news)} new rows, {missing} without a published row, {len(cells)} cells')
print(f"  cells with all 500 replicates equal: {sum(x['reps'] == 500 and x['equal'] == 500 for x in v3)} of {len(v3)}")
for x in v3:
    if not (x['reps'] == 500 and x['equal'] == 500):
        print('  NOT EQUAL:', x)
# medians by n over the 12 cells (2 dists x 2 d x 3 k)
by_n = defaultdict(lambda: {'old': [], 'new': []})
for x in v3:
    by_n[x['n']]['old'].append(x['median_sec_old']); by_n[x['n']]['new'].append(x['median_sec_new'])
summ = {}
for n in sorted(by_n):
    o, w = by_n[n]['old'], by_n[n]['new']
    summ[n] = {'cells': len(o), 'old_min': min(o), 'old_max': max(o), 'old_median': statistics.median(o),
               'new_min': min(w), 'new_max': max(w), 'new_median': statistics.median(w)}
    print(f"  n={n:3d}: median seconds per data set over {len(o)} cells: old {min(o):.3f}-{max(o):.3f} "
          f"(median {statistics.median(o):.3f}); new {min(w):.3f}-{max(w):.3f} (median {statistics.median(w):.3f})")
if 200 in summ and 500 in summ:
    for t in ('old', 'new'):
        sl = math.log(summ[500][f'{t}_median'] / summ[200][f'{t}_median']) / math.log(500 / 200)
        summ[f'slope_200_500_{t}'] = sl
        print(f'  log-log slope n=200 -> 500 ({t}): {sl:.2f}')
out['v3_by_n'] = summ
# as the paper's S3 numbers were computed: per (dist, d), the median over the 1500 data sets (3 k x 500)
# at n = 200 and n = 500, and the exponent log(m500 / m200) / log(2.5) per (dist, d)
pool_o = defaultdict(list); pool_n = defaultdict(list)
for c in sec_old:
    pool_o[(c[0], c[1], c[2])] += sec_old[c]; pool_n[(c[0], c[1], c[2])] += sec_new[c]
grp = {}
print('\n  per (dist, d), pooled over k: median seconds n=200 | n=500 | exponent')
for dist in ('uniform', 'gaussian'):
    for d in ('10', '20'):
        g = {}
        for t, pool in (('old', pool_o), ('new', pool_n)):
            if (dist, d, '200') in pool and (dist, d, '500') in pool:
                a = statistics.median(pool[(dist, d, '200')]); b = statistics.median(pool[(dist, d, '500')])
                g[t] = {'n200': a, 'n500': b, 'exp': math.log(b / a) / math.log(2.5),
                        'count200': len(pool[(dist, d, '200')]), 'count500': len(pool[(dist, d, '500')])}
        grp[f'{dist}_d{d}'] = g
        print(f"  {dist:8s} d={d}: " + '; '.join(f"{t} {v['n200']:.3f} | {v['n500']:.3f} | {v['exp']:.2f} "
                                                   f"(n={v['count200']}/{v['count500']})" for t, v in g.items()))
for t in ('old', 'new'):
    vals = [g[t] for g in grp.values() if t in g]
    if vals:
        print(f"  range {t}: n=200 {min(v['n200'] for v in vals):.2f}-{max(v['n200'] for v in vals):.2f} s, "
              f"n=500 {min(v['n500'] for v in vals):.2f}-{max(v['n500'] for v in vals):.2f} s, "
              f"exponent {min(v['exp'] for v in vals):.2f}-{max(v['exp'] for v in vals):.2f}")
out['v3_pooled'] = grp
json.dump(out, open(os.path.join(S, 'compare_summary.json'), 'w'), indent=1)
print('\nwritten', os.path.join(S, 'compare_summary.json'))
