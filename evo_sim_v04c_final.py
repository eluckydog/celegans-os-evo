#!/usr/bin/env python3
"""
C. elegans OS 演化模拟器 v0.4c — 30参数 + 白盒 + 跨层约束
"""
import random, json, csv, os, time
import numpy as np
from pathlib import Path

ROOT = Path(__file__).parent
RESULTS_DIR = ROOT / "results_v04"
RESULTS_DIR.mkdir(exist_ok=True)

WILD_PARAMS = {
    "lin12_expr": 1.0, "lag2_expr": 1.0, "notch_sens": 1.0,
    "notch_thresh": 0.5, "lat_inhib": 1.0, "notch_decay": 0.5, "dsl_basal": 0.5,
    "fork_rate": 1.0, "fork_accuracy": 0.99, "mother_size": 0.8, "daughter_asym": 0.5, "polarity": 0.8,
    "timer_speed": 1.0, "timer_precis": 0.8, "check_g1": 0.7, "check_g2": 0.7, "dna_damage": 0.8,
    "ach_sens": 0.6, "gaba_sens": 0.6, "glu_sens": 0.5, "serotonin": 0.5,
    "memory_stab": 0.9, "h3k4me3": 0.7, "h3k27me3": 0.7, "chrom_noise": 0.3, "heritability": 0.8,
    "apop_thresh": 0.6, "caspase": 0.5, "egl1_sens": 0.5, "cell_survive": 0.8,
}

COUPLING_RULES = [
    ("lag2_expr", "lin12_expr", 0.05), ("notch_decay", "notch_sens", 0.03),
    ("fork_rate", "timer_speed", 0.04), ("cell_survive", "apop_thresh", 0.03),
    ("caspase", "egl1_sens", 0.03), ("h3k4me3", "h3k27me3", -0.03),
    ("memory_stab", "heritability", 0.03), ("check_g1", "check_g2", 0.03),
    ("mother_size", "daughter_asym", 0.02), ("ach_sens", "gaba_sens", 0.02),
]


class GeneNetwork:
    def __init__(self, params=None):
        self.params = {k: v for k, v in WILD_PARAMS.items()}
        if params: self.params.update(params)

    def mutate(self, rate=0.12, severity=0.15):
        c = GeneNetwork()
        c.params = {k: v for k, v in self.params.items()}
        for k in c.params:
            if random.random() < rate:
                p = random.random()
                if p < 0.8: c.params[k] *= (1 + np.random.normal(0, severity))
                elif p < 0.9: c.params[k] *= 0.1
                else: c.params[k] *= 3.0
                c.params[k] = max(0.01, min(10.0, c.params[k]))
        return c

    def crossover(self, other):
        c = GeneNetwork()
        for k in c.params: c.params[k] = random.choice([self.params[k], other.params[k]])
        return c


def coupling_penalty(params):
    penalty = 0.0
    for k, v in params.items():
        r = v / WILD_PARAMS[k]
        if r > 5.0: penalty += (r - 5.0) * 0.3
        if r < 0.1: penalty += (0.1 - r) * 0.5
    for a, b, w in COUPLING_RULES:
        va = params[a] / WILD_PARAMS[a] - 1.0
        vb = params[b] / WILD_PARAMS[b] - 1.0
        expectation = w * va * vb
        penalty -= max(-0.1, min(0.1, expectation))
    return penalty


# --- 子任务 ---

def task_ac_vu(params, n_seeds=3):
    scores = []
    for seed in range(n_seeds):
        rng = np.random.RandomState(seed)
        N1, N2 = 0.1+rng.rand()*0.01, 0.1+rng.rand()*0.01
        G1, G2 = 0.5, 0.5
        alpha, beta, gamma = 2.0, 1.0, 0.2
        K_base, dt = 5.0, 0.1
        N1_p, N2_p = N1, N2
        for step in range(250):
            raw1 = G2 * params["notch_sens"] * params["dsl_basal"]
            raw2 = G1 * params["notch_sens"] * params["dsl_basal"]
            N1 += dt * (params["lin12_expr"]*alpha*raw1*(1-N1) - params["notch_decay"]*N1 + 0.001*(rng.rand()-0.5))
            N2 += dt * (params["lin12_expr"]*alpha*raw2*(1-N2) - params["notch_decay"]*N2)
            N1, N2 = max(0,min(2,N1)), max(0,min(2,N2))
            inh1 = 1.0 + K_base*N1**2/(N1**2+0.1)
            inh2 = 1.0 + K_base*N2**2/(N2**2+0.1)
            G1 += dt * (params["lag2_expr"]*beta/inh1 - gamma*G1)
            G2 += dt * (params["lag2_expr"]*beta/inh2 - gamma*G2)
            G1, G2 = max(0,min(5,G1)), max(0,min(5,G2))
            if step > 30 and abs(N1-N2) > 0.3:
                if abs(N1-N1_p) < 0.001 and abs(N2-N2_p) < 0.001: break
            N1_p, N2_p = N1, N2
        diff = abs(N1-N2)
        s = 1.0 + min(0.3, diff*0.5) if diff > 0.05 else 0.0
        scores.append(s)
    return np.mean(scores)


def task_fork(params, n_seeds=3):
    scores = []
    for seed in range(n_seeds):
        rng = np.random.RandomState(seed)
        n_div, clock, t = 0, 0.0, 0.0
        br = params["fork_rate"] * (1.0 + params["polarity"]*0.2)
        while n_div < 5 and t < 500:
            clock += br * 0.2; t += 1
            if clock >= 1.0:
                clock = 0.0; n_div += 1
                if rng.rand() > params["fork_accuracy"]: n_div += 0.3-rng.rand()
        s = min(1.0, n_div/5.0)
        if t < 100: s += 0.2
        scores.append(s)
    return np.mean(scores)


def task_timer(params, n_seeds=3):
    scores = []
    for seed in range(n_seeds):
        rng = np.random.RandomState(seed)
        cycles = []
        for _ in range(10):
            t = params["timer_speed"]
            if rng.rand() < (1-params["timer_precis"]): t *= 0.5
            if rng.rand() < (1-params["check_g1"]): t *= 0.7
            if rng.rand() < (1-params["check_g2"]): t *= 0.7
            if rng.rand() < (1-params["dna_damage"]): t *= 0.5
            cycles.append(t)
        avg, sd = np.mean(cycles), np.std(cycles)
        s = min(1.0, avg*0.5)
        if sd < 0.2: s += 0.2
        if avg > 0.5: s += 0.1
        scores.append(s)
    return np.mean(scores)


def task_noise(params, n_seeds=10):
    stab = params["ach_sens"]*0.3 + params["gaba_sens"]*0.4 + params["glu_sens"]*0.2 + params["serotonin"]*0.1
    rng = np.random.RandomState(42)
    sc = []
    for lv in [0.001, 0.005, 0.01, 0.05, 0.1]:
        N = 0.5 + rng.rand()*0.1
        for _ in range(100):
            N += lv*(rng.rand()-0.5)/max(0.1,stab) - N*0.1
            N = max(0,min(1,N))
        sc.append(max(0,1-abs(N-0.5)*2))
    return np.mean(sc)


def task_memory(params, n_seeds=5):
    scores = []
    for seed in range(n_seeds):
        rng = np.random.RandomState(seed)
        state = 0.0
        decay = max(0.01, 0.5-params["memory_stab"]*0.4)
        rate = params["h3k4me3"]*0.3 + params["memory_stab"]*0.3
        repress = params["h3k27me3"]*0.3
        for step in range(100):
            if step < 10: state += rate*0.1
            state -= state*decay
            state += rng.randn()*params["chrom_noise"]*0.1
            state -= repress*0.02
            state = max(0,min(1,state))
        scores.append(max(0, min(1, state*2 + params["heritability"]*0.3)))
    return np.mean(scores)


def task_apoptosis(params, n_seeds=5):
    scores = []
    for stress in [0.1, 0.3, 0.5, 0.7, 0.9]:
        dp = stress * (1.0 - params["cell_survive"]*0.5)
        if dp > params["apop_thresh"]: dp *= (1 + params["caspase"]*0.3)
        if params["egl1_sens"] > 0.5: dp *= (1 + (params["egl1_sens"]-0.5)*0.5)
        scores.append(1.0 - dp)
    return np.mean(scores)


TASKS = {"ac_vu": task_ac_vu, "fork": task_fork, "timer": task_timer,
         "noise": task_noise, "memory": task_memory, "apoptosis": task_apoptosis}
WEIGHTS = {"ac_vu": 1.0, "fork": 0.6, "timer": 0.5, "noise": 0.4, "memory": 0.4, "apoptosis": 0.4}


def fitness_30(net, n_seeds=3, verbose=False):
    total = 0.0
    sc = {}
    for name, fn in TASKS.items():
        s = fn(net.params, n_seeds=n_seeds)
        sc[name] = round(s, 3)
        total += s * WEIGHTS[name]
    total -= coupling_penalty(net.params) * 0.5
    return max(0.0, total)


def create_pop(pop_size=50):
    pop = []
    for _ in range(pop_size):
        ind = GeneNetwork()
        for k in ind.params: ind.params[k] *= (1+np.random.normal(0,0.1))
        pop.append(ind)
    return pop


def detect_tips(traj):
    pts = []
    for i in range(1, len(traj)):
        db = traj[i]["best"] - traj[i-1]["best"]
        da = traj[i]["avg"] - traj[i-1]["avg"]
        if db >= 0.15 or db+da >= 0.2:
            pts.append({"gen": traj[i]["gen"], "db": round(db,3), "da": round(da,3)})
    return pts


def evolve(label, n_gens=100, pop_size=50, mut_rate=0.12, n_seeds=3):
    pop = create_pop(pop_size)
    traj = []
    start = time.time()

    for gen in range(n_gens):
        scores = [fitness_30(ind, n_seeds=n_seeds) for ind in pop]
        bi = np.argmax(scores)
        best_diff = {k: round(pop[bi].params[k]-WILD_PARAMS[k], 3)
                    for k in WILD_PARAMS if abs(pop[bi].params[k]-WILD_PARAMS[k]) > 0.1}
        mean_p = {k: round(float(np.mean([ind.params[k] for ind in pop])), 3) for k in ["lin12_expr","lag2_expr","fork_rate","cell_survive","egl1_sens","notch_thresh"]}

        if gen % 10 == 0 or gen == n_gens-1:
            bp = "; ".join(f"{k}={v}" for k,v in sorted(best_diff.items())[:5])
            mp = "; ".join(f"{k}={v}" for k,v in mean_p.items())
            print(f"  Gen {gen:3d} | b={scores[bi]:.3f} a={np.mean(scores):.3f} | pop_mean: {mp}")

        traj.append({
            "gen": gen, "best": scores[bi], "avg": float(np.mean(scores)),
            "worst": float(np.min(scores)), "std": float(np.std(scores)),
            "best_params": best_diff, "mean_params": mean_p,
        })

        ranked = sorted(range(len(scores)), key=lambda i: -scores[i])
        ranked_inds = [(scores[i], pop[i]) for i in ranked]
        n_elite = max(1, int(pop_size*0.1))
        new_pop = [ranked_inds[i][1] for i in range(n_elite)]
        while len(new_pop) < pop_size:
            pool = random.choices(ranked_inds, k=3)
            p1 = max(pool[:2], key=lambda x: x[0])[1]
            p2 = max(pool[1:3], key=lambda x: x[0])[1]
            child = p1.crossover(p2).mutate(mut_rate, 0.15)
            new_pop.append(child)
        pop = new_pop

    elapsed = time.time() - start
    tips = detect_tips(traj)

    print(f"\n  Done: {elapsed:.1f}s | b={traj[-1]['best']:.3f} a={traj[-1]['avg']:.3f}")
    print(f"  Tips: {len(tips)} points")

    return pop, traj, tips, elapsed


def save_outputs(pop, traj, tips, elapsed, label):
    csv_path = RESULTS_DIR / f"v04c_{label}.csv"
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["gen","best","avg","worst","std","lin12","lag2","fork_rate","cell_survive","egl1_sens","notch_thresh"])
        for s in traj:
            mp = s.get("mean_params", {})
            w.writerow([s["gen"],s["best"],s["avg"],s["worst"],s["std"],
                       mp.get("lin12_expr",""), mp.get("lag2_expr",""),
                       mp.get("fork_rate",""), mp.get("cell_survive",""),
                       mp.get("egl1_sens",""), mp.get("notch_thresh","")])

    rpt = {
        "label": label, "time": round(elapsed,1), "gens": len(traj),
        "best": traj[-1]["best"], "avg": traj[-1]["avg"],
        "tips": len(tips),
        "best_diff": {k: round(pop[0].params[k]-WILD_PARAMS[k],3)
                     for k in WILD_PARAMS if abs(pop[0].params[k]-WILD_PARAMS[k]) > 0.15},
    }
    rpt_path = RESULTS_DIR / f"v04c_{label}_report.json"
    with open(rpt_path, 'w') as f:
        json.dump(rpt, f, indent=2)
    return csv_path, rpt_path


if __name__ == "__main__":
    print(f"{'='*65}")
    print(f"  C. elegans OS v0.4c — 30参数 + 白盒 + 跨层约束")
    print(f"  {len(COUPLING_RULES)}耦合规则 | 6子任务 | 4方向 x 100代")
    print(f"{'='*65}")

    dirs = {
        "baseline": {"mut_rate": 0.12, "pop_size": 50},
        "low_mut": {"mut_rate": 0.05, "pop_size": 50},
        "high_mut": {"mut_rate": 0.20, "pop_size": 50},
        "small_pop": {"mut_rate": 0.12, "pop_size": 20},
    }

    for key, cfg in dirs.items():
        print(f"\n{'='*60}")
        print(f"  Dir: {key} ({cfg['pop_size']}pop, mut={cfg['mut_rate']})")
        print(f"{'='*60}")
        pop, traj, tips, elapsed = evolve(key, n_gens=100, **cfg)
        csv_path, rpt_path = save_outputs(pop, traj, tips, elapsed, key)
        print(f"  CSV: {csv_path.name} | RPT: {rpt_path.name}")

    print(f"\n{'='*65}")
    print(f"  All done. Results in: {RESULTS_DIR}")
    print(f"{'='*65}")
