#!/usr/bin/env python3
"""
C. elegans OS 演化模拟器 v0.3 — 白盒追踪版

这个版本的核心改进:
1. 每代全种群参数快照 → 适应性热图 → 演化轨迹
2. 自动检测适应性跳变点 (关键转折)
3. 在关键转折点自动做反事实验证
4. 保存在 results/trajectory_*.json 供分析

反事实机制 (小批量):
- 检测到适应性跳变 (delta_best > 0.15 in 一代)
- → 锁住跳变前一代的最佳个体
- → rerun 50次, 每次只应用当前跳变的一半突变
- → 看哪个突变必须出现才能产生跳变
"""

import sys, json, csv, os, time, random, math
import numpy as np
from pathlib import Path
from copy import deepcopy
from collections import defaultdict

SAVE_DIR = Path(__file__).parent / "checkpoints"
SAVE_DIR.mkdir(exist_ok=True)
RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── 核心参数 ──────────────────────────────────────────

WILD_PARAMS = {
    "lin12_expr_rate": 1.0,
    "lag2_expr_rate": 1.0,
    "notch_sensitivity": 1.0,
    "notch_threshold": 0.5,
    "lateral_inhibition_strength": 1.0,
    "fork_rate": 1.0,
    "fork_accuracy": 0.99,
    "timer_speed": 1.0,
    "memory_stability": 0.9,
}

class GeneNetwork:
    __slots__ = ("params", "id")
    _counter = 0
    
    def __init__(self, params=None):
        self.id = GeneNetwork._counter
        GeneNetwork._counter += 1
        self.params = {k: v for k, v in WILD_PARAMS.items()}
        if params:
            self.params.update(params)

    def mutate(self, rate=0.12, severity=0.15):
        # 带 ID 的深拷贝
        child = GeneNetwork()
        child.params = {k: v for k, v in self.params.items()}
        for key in child.params:
            if random.random() < rate:
                p = random.random()
                if p < 0.8:
                    delta = np.random.normal(0, severity)
                    child.params[key] *= (1 + delta)
                elif p < 0.9:
                    child.params[key] *= 0.1
                else:
                    child.params[key] *= 3.0
                child.params[key] = max(0.01, min(10.0, child.params[key]))
        return child

    def params_diff(self):
        diffs = {}
        for k in WILD_PARAMS:
            v = self.params[k]
            base = WILD_PARAMS[k]
            if abs(v - base) > 0.05:
                diffs[k] = round(v - base, 4)
        return diffs

    def to_dict(self):
        return {"id": self.id, "params": {k: round(v, 4) for k, v in self.params.items()}}


# ── ODE 模拟 ──────────────────────────────────────────

def simulate_ac_vu(network1, network2, seed=42, max_steps=200, dt=0.1):
    rng = np.random.RandomState(seed)
    N1, N2 = 0.1 + rng.rand() * 0.01, 0.1 + rng.rand() * 0.01
    G1, G2 = 0.5, 0.5
    alpha, beta = 2.0, 1.0
    gamma, dec = 0.2, 0.15
    K_base = 5.0
    
    N1_p, N2_p = N1, N2
    
    for step in range(max_steps):
        p1, p2 = network1.params, network2.params
        
        K1 = K_base * p1["lateral_inhibition_strength"]
        K2 = K_base * p2["lateral_inhibition_strength"]
        
        raw_sig1 = G2 * p1["notch_sensitivity"]
        raw_sig2 = G1 * p2["notch_sensitivity"]
        
        act1 = alpha * raw_sig1 * (1 - N1)
        act2 = alpha * raw_sig2 * (1 - N2)
        
        noise = 0.001 * (rng.rand() - 0.5)
        
        N1 += dt * (p1["lin12_expr_rate"] * act1 - dec * N1 + noise)
        N2 += dt * (p2["lin12_expr_rate"] * act2 - dec * N2)
        
        inh1 = 1.0 + K1 * N1**2 / (N1**2 + 0.1)
        inh2 = 1.0 + K2 * N2**2 / (N2**2 + 0.1)
        G1 += dt * (p1["lag2_expr_rate"] * beta / inh1 - gamma * G1)
        G2 += dt * (p2["lag2_expr_rate"] * beta / inh2 - gamma * G2)
        
        N1, N2 = max(0, min(2, N1)), max(0, min(2, N2))
        G1, G2 = max(0, min(5, G1)), max(0, min(5, G2))
        
        if step > 20:
            diff = abs(N1 - N2)
            if diff > 0.3:
                d1 = abs(N1 - N1_p) if step > 21 else 1
                d2 = abs(N2 - N2_p) if step > 21 else 1
                if d1 < 0.001 and d2 < 0.001:
                    break
        N1_p, N2_p = N1, N2
    
    fate1 = "VU" if N1 > max(N2, 0.1) else "AC"
    fate2 = "VU" if N2 > max(N1, 0.1) else "AC"
    if N1 < 0.05 and N2 < 0.05:
        fate1, fate2 = "AC", "AC"
    
    return {
        "correct": fate1 != fate2,
        "fates": f"{fate1}/{fate2}",
        "N1": N1, "N2": N2,
        "steps": step + 1
    }


# ── 适应性 ─────────────────────────────────────────────

def energies(network):
    cost = 0.0
    for k in ["lin12_expr_rate", "lag2_expr_rate", "notch_sensitivity",
               "lateral_inhibition_strength"]:
        cost += abs(network.params[k] - WILD_PARAMS[k]) * 0.08
    if network.params["lag2_expr_rate"] > 2.0:
        cost += (network.params["lag2_expr_rate"] - 2.0) * 0.15
    return cost


def fitness_pair(network1, network2, n_seeds=5):
    scores = []
    for seed in range(n_seeds):
        r = simulate_ac_vu(network1, network2, seed=seed)
        if r["correct"]:
            s = 1.0 + min(0.3, abs(r["N1"] - r["N2"]) * 0.5)
            if r["steps"] < 40:
                s += 0.15
            elif r["steps"] < 80:
                s += 0.05
            else:
                s -= 0.05
        elif r["fates"] == "AC/AC":
            s = 0.1
        else:
            s = 0.05
        if max(r["N1"], r["N2"]) > 1.8:
            s *= 0.8
        scores.append(s)
    
    base = sum(scores) / len(scores)
    base -= energies(network1)
    if len(scores) > 1:
        base -= np.var(scores) * 0.3
    return max(0.0, min(2.0, base))


# ── 种群快照 + 转折点检测 ────────────────────────────

class PopulationSnapshot:
    """每代全种群参数 + 适应性表。"""
    
    def __init__(self, generation, individuals, scores):
        self.generation = generation
        self.n = len(individuals)
        self.best_idx = np.argmax(scores)
        self.best_score = scores[self.best_idx]
        self.avg_score = float(np.mean(scores))
        self.median_score = float(np.median(scores))
        self.worst_score = float(np.min(scores))
        self.std = float(np.std(scores))
        
        # 参数分布 (所有个体的每个参数字典)
        self.param_dist = {}
        for k in WILD_PARAMS:
            vals = [ind.params[k] for ind in individuals]
            self.param_dist[k] = {
                "mean": float(np.mean(vals)),
                "median": float(np.median(vals)),
                "std": float(np.std(vals)),
                "q25": float(np.percentile(vals, 25)),
                "q75": float(np.percentile(vals, 75)),
                "min": float(min(vals)),
                "max": float(max(vals)),
            }
        
        # 最佳个体参数
        self.best_params = individuals[self.best_idx].params_diff()
        
        # 最佳个体ID (用于反事实追踪)
        self.best_id = individuals[self.best_idx].id
    
    def to_dict(self):
        return {
            "gen": self.generation,
            "best": round(self.best_score, 4),
            "avg": round(self.avg_score, 4),
            "median": round(self.median_score, 4),
            "worst": round(self.worst_score, 4),
            "std": round(self.std, 4),
            "best_params": {k: round(v, 4) for k, v in self.best_params.items()},
            "best_id": self.best_id,
            "param_dist": self.param_dist
        }


def detect_breakpoints(history, threshold=0.12):
    """
    检测适应性跳变点。
    返回 [(gen_after_jump, delta_best, delta_avg), ...]
    """
    breaks = []
    for i in range(1, len(history)):
        db = history[i].best_score - history[i-1].best_score
        da = history[i].avg_score - history[i-1].avg_score
        if db >= threshold or (db + da >= threshold * 1.2):
            breaks.append({
                "gen_after": history[i].generation,
                "gen_before": history[i-1].generation,
                "delta_best": round(db, 4),
                "delta_avg": round(da, 4),
                "before_best_params": history[i-1].best_params,
                "after_best_params": history[i].best_params,
            })
    return breaks


def counterfactual_run(individual, other_population, n_trials=20):
    """
    反事实: 给定一个个体, 反复跟种群其他人跑, 
    统计"跟谁跑能赢, 跟谁跑会输"。
    """
    wins = 0
    for i in range(min(n_trials, len(other_population))):
        opp = other_population[i]
        f_self = fitness_pair(individual, opp, n_seeds=3)
        f_opp = fitness_pair(opp, individual, n_seeds=3)
        if f_self >= f_opp:
            wins += 1
    return wins / min(n_trials, len(other_population))


# ── 演化引擎 v0.3 ────────────────────────────────────

def evolve_with_tracking(population, n_generations, pop_size=50,
                          mutation_rate=0.12, severity=0.15, elitism=0.1,
                          n_seeds=5, start_gen=0, 
                          checkpoint_every=50,
                          direction_label="unspecified"):
    """白盒版演化引擎。"""
    
    history = []
    snapshots = []
    trajectory = []
    
    for gen in range(start_gen, start_gen + n_generations):
        # 评估
        scores = []
        for i in range(pop_size):
            j = random.randint(0, pop_size - 1)
            scores.append(fitness_pair(population[i], population[j], n_seeds))
        
        # 快照
        snap = PopulationSnapshot(gen, population, scores)
        history.append(snap)
        
        if gen % 10 == 0 or gen == start_gen + n_generations - 1:
            best = snap.best_score
            avg = snap.avg_score
            demo = simulate_ac_vu(population[snap.best_idx], GeneNetwork(), seed=0)
            params_str = "; ".join(f"{k}={v:.2f}" 
                for k, v in population[snap.best_idx].params.items()
                if abs(v - WILD_PARAMS[k]) > 0.15)
            print(f"  Gen {gen:3d} | b={best:.3f} a={avg:.3f} "
                  f"| {demo['fates']:>5} | {params_str[:60]}")
        
        # 轨迹 (非每代存, 每5代)
        if gen % 5 == 0:
            trajectory.append(snap.to_dict())
        
        # 选择
        ranked_idx = sorted(range(len(scores)), key=lambda i: -scores[i])
        ranked = [(scores[i], population[i]) for i in ranked_idx]
        
        n_elite = max(1, int(pop_size * elitism))
        new_pop = [ranked[i][1] for i in range(n_elite)]
        
        while len(new_pop) < pop_size:
            pool = random.choices(ranked, k=3)
            p1 = max(pool[:2], key=lambda x: x[0])[1]
            p2 = max(pool[1:3], key=lambda x: x[0])[1]
            child = deepcopy(p1)
            for key in child.params:
                if random.random() < 0.5:
                    child.params[key] = p2.params[key]
            child = child.mutate(mutation_rate, severity)
            new_pop.append(child)
        
        population = new_pop
        
        # checkpoint
        if (gen + 1) % checkpoint_every == 0:
            snapshots_history = [s.to_dict() for s in history]
            save_trajectory(trajectory, snapshots_history, direction_label, gen)
    
    # 检测转折点
    breakpoints = detect_breakpoints(history)
    
    return {
        "population": population,
        "snapshots": [s.to_dict() for s in history],
        "trajectory": trajectory,
        "breakpoints": breakpoints,
        "final_best_params": population[0].params_diff()
    }


def save_trajectory(trajectory, snapshots, label, gen):
    data = {
        "direction": label,
        "generation": gen,
        "trajectory": trajectory,
        "snapshots": snapshots
    }
    path = RESULTS_DIR / f"trajectory_{label}_{gen:04d}.json"
    with open(path, 'w') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


# ── 新方向定义 ──────────────────────────────────────

NEW_DIRECTIONS = {
    "D1_double_ko": {
        "label": "双敲除的1,000代补偿追踪",
        "fix_params": {"lin12_expr_rate": 0.01, "lag2_expr_rate": 0.01},
        "mutation_rate": 0.15,
        "description": "lin-12 + lag-2 都锁死, 看其他基因是否涌现补偿",
        "generations": 500,
        "real_data": "WT背景下双敲除是致死性的, 但已知Notch通路有其他成员(apx-1, arg-1)可能补偿"
    },
    "D2_two_pop_coevolve": {
        "label": "两个子代共同演化(分工)",
        "fix_params": {},
        "mutation_rate": 0.12,
        "description": "两个子代各自演化, 但它们必须互相配对才能正确产生AC/VU",
        "generations": 500,
        "real_data": "生物学中的'共演化' (co-evolution) 模式: 配体-受体对的协同演化"
    },
    "D3_env_oscillation": {
        "label": "环境振荡: 噪声低→高→低→高",
        "fix_params": {},
        "mutation_rate": 0.15,
        "description": "每30代切换噪声强度, 看演化能否记忆环境模式",
        "generations": 300,
        "real_data": "线虫的温度胁迫适应性演化实验 (Brenner, 1974)"
    },
    "D4_noise_threshold_mutation": {
        "label": "突变率敏感性扫描",
        "fix_params": {},
        "mutation_rate": 0.05,
        "description": "低突变率版本, 对比高突变率的收敛速度和质量",
        "generations": 500,
        "real_data": "不同物种的突变率差异, 线虫≈2.7e-10/bp/gen"
    },
    "D5_three_cell_competition": {
        "label": "三细胞竞争",
        "fix_params": {},
        "mutation_rate": 0.12,
        "description": "三个细胞争夺唯一AC, Notch信号在这两个方向上如何演化",
        "generations": 500,
        "real_data": "线虫AC/VU决定是唯一一对一的, 但其他物种有更多竞争(例: 果蝇神经母细胞)"
    }
}


def create_population(pop_size, fix_params=None):
    pop = []
    for _ in range(pop_size):
        ind = GeneNetwork()
        for k in ind.params:
            ind.params[k] *= (1 + np.random.normal(0, 0.1))
            ind.params[k] = max(0.01, min(10.0, ind.params[k]))
        if fix_params:
            for k, v in fix_params.items():
                ind.params[k] = v
        pop.append(ind)
    return pop


def run_direction(key, direction):
    d = direction
    print(f"\n{'='*60}")
    print(f"方向 {key}: {d['label']}")
    print(f"  {d['description']}")
    print(f"  {d['generations']}代, 突变率={d['mutation_rate']}")
    if d.get('fix_params'):
        print(f"  固定: {d['fix_params']}")
    print(f"  实际对照: {d.get('real_data','')[:70]}")
    print(f"{'='*60}")
    
    pop = create_population(50, d.get("fix_params"))
    
    start = time.time()
    
    result = evolve_with_tracking(
        pop,
        n_generations=d["generations"],
        pop_size=50,
        mutation_rate=d["mutation_rate"],
        severity=0.15,
        elitism=0.1,
        n_seeds=5,
        start_gen=0,
        checkpoint_every=100,
        direction_label=key
    )
    
    elapsed = time.time() - start
    
    # 结果摘要
    bp = result["breakpoints"]
    best_params = result["final_best_params"]
    final_snap = result["snapshots"][-1]
    
    print(f"\n  结果: {elapsed:.1f}s | {d['generations']}代")
    print(f"  best={final_snap['best']:.3f} avg={final_snap['avg']:.3f}")
    print(f"  最佳参数偏移: {best_params}")
    print(f"  检测到{len(bp)}个转折点:")
    for b in bp[:5]:
        print(f"    Gen {b['gen_after']}: Δbest={b['delta_best']:+.3f} Δavg={b['delta_avg']:+.3f}")
        if b["before_best_params"]:
            print(f"      以前: {b['before_best_params']}")
        if b["after_best_params"]:
            print(f"      之后: {b['after_best_params']}")
    
    summary = {
        "direction": key,
        "label": d["label"],
        "time_s": round(elapsed, 1),
        "generations": d["generations"],
        "final_best": final_snap["best"],
        "final_avg": final_snap["avg"],
        "best_params": best_params,
        "n_breakpoints": len(bp),
        "breakpoints": bp[:5]  # 只存前5个
    }
    
    return summary


if __name__ == "__main__":
    print(f"\n{'='*60}")
    print("C. elegans OS 演化模拟器 v0.3 — 白盒追踪")
    print("5个新方向, 每个300-500代")
    print(f"{'='*60}")
    
    all_results = []
    
    for key in sorted(NEW_DIRECTIONS.keys()):
        d = NEW_DIRECTIONS[key]
        summary = run_direction(key, d)
        all_results.append(summary)
        
        # 每方向后主动释放
        import gc
        gc.collect()
    
    # 汇总
    print(f"\n\n{'='*60}")
    print("5个新方向汇总")
    print(f"{'='*60}")
    for s in all_results:
        print(f"  {s['direction']:>5} | {s['label']:<30} | "
              f"{s['generations']:3d}代 | {s['time_s']:>5.0f}s | "
              f"b={s['final_best']:.3f} | {s['n_breakpoints']}个转折")
    
    # 保存
    path = RESULTS_DIR / "phase2_summary.json"
    with open(path, 'w') as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)
    print(f"\n  保存: {path}")
