#!/usr/bin/env python3
"""
C. elegans OS 演化模拟器 v0.4 — 30参数 + 多任务适应性

覆盖7个OS层:
- Layer 2 (IPC-Notch)   — 7个参数
- Layer 3 (fork)        — 5个参数
- Layer 4 (调度器)      — 5个参数
- Layer 5 (驱动-神经递质) — 4个参数
- Layer 6 (染色质)      — 5个参数
- Layer 7 (凋亡)        — 4个参数

多任务适应性:
  1. AC/VU正确性 (0-1.0分, 核心)
  2. 发育速度 (0-0.3分)
  3. 噪声容忍 (0-0.2分)  
  4. 凋亡正确性 (0-0.2分)
  5. 能量成本 (扣分项)
"""

import random, json, csv, os, time
import numpy as np
from pathlib import Path
from copy import deepcopy
from collections import defaultdict

# ── 文件路径 ────────────────────────────────────────────

ROOT = Path(__file__).parent
SAVE_DIR = ROOT / "checkpoints_v04"
SAVE_DIR.mkdir(exist_ok=True)
RESULTS_DIR = ROOT / "results_v04"
RESULTS_DIR.mkdir(exist_ok=True)

# ── 30参数定义 ──────────────────────────────────────────

WILD_PARAMS = {
    # Layer 2: IPC-Notch (7)
    "lin12_expr": 1.0,
    "lag2_expr": 1.0,
    "notch_sens": 1.0,
    "notch_thresh": 0.5,
    "lat_inhib": 1.0,
    "notch_decay": 0.5,
    "dsl_basal": 0.5,
    # Layer 3: fork (5)
    "fork_rate": 1.0,
    "fork_accuracy": 0.99,
    "mother_size": 0.8,
    "daughter_asym": 0.5,
    "polarity": 0.8,
    # Layer 4: scheduler/timer (5)
    "timer_speed": 1.0,
    "timer_precis": 0.8,
    "check_g1": 0.7,
    "check_g2": 0.7,
    "dna_damage": 0.8,
    # Layer 5: drivers (4)
    "ach_sens": 0.6,
    "gaba_sens": 0.6,
    "glu_sens": 0.5,
    "serotonin": 0.5,
    # Layer 6: chromatin (5)
    "memory_stab": 0.9,
    "h3k4me3": 0.7,
    "h3k27me3": 0.7,
    "chrom_noise": 0.3,
    "heritability": 0.8,
    # Layer 7: apoptosis (4)
    "apop_thresh": 0.6,
    "caspase": 0.5,
    "egl1_sens": 0.5,
    "cell_survive": 0.8,
}

PARAM_GROUPS = {
    "Notch": ["lin12_expr", "lag2_expr", "notch_sens", "notch_thresh", "lat_inhib", "notch_decay", "dsl_basal"],
    "Fork": ["fork_rate", "fork_accuracy", "mother_size", "daughter_asym", "polarity"],
    "Timer": ["timer_speed", "timer_precis", "check_g1", "check_g2", "dna_damage"],
    "Drivers": ["ach_sens", "gaba_sens", "glu_sens", "serotonin"],
    "Chromatin": ["memory_stab", "h3k4me3", "h3k27me3", "chrom_noise", "heritability"],
    "Apoptosis": ["apop_thresh", "caspase", "egl1_sens", "cell_survive"],
}


class GeneNetwork:
    __slots__ = ("params",)
    
    def __init__(self, params=None):
        self.params = {k: v for k, v in WILD_PARAMS.items()}
        if params:
            self.params.update(params)
    
    def mutate(self, rate=0.12, severity=0.15):
        child = GeneNetwork()
        child.params = {k: v for k, v in self.params.items()}
        for key in child.params:
            if random.random() < rate:
                p = random.random()
                if p < 0.8:
                    child.params[key] *= (1 + np.random.normal(0, severity))
                elif p < 0.9:
                    child.params[key] *= 0.1
                else:
                    child.params[key] *= 3.0
                child.params[key] = max(0.01, min(10.0, child.params[key]))
        return child
    
    def crossover(self, other):
        child = GeneNetwork()
        for key in child.params:
            child.params[key] = random.choice([self.params[key], other.params[key]])
        return child


# ── ODE: 3细胞Notch + 时钟 + 凋亡 ──────────────────────

def simulate_multi_cell(network, n_cells=3, seed=42, max_steps=500, dt=0.05):
    """
    多细胞Notch竞争模拟。
    n_cells个细胞初始对称, 通过Notch侧向抑制竞争一个"AC"命运。
    同时运行: 发育时钟, 凋亡信号。
    """
    rng = np.random.RandomState(seed)
    
    p = network.params
    n = n_cells
    
    # Notch信号变量
    N = np.array([0.1 + rng.rand() * 0.01 for _ in range(n)])
    G = np.array([0.5 for _ in range(n)])
    
    # 发育时钟
    cell_cycle = np.zeros(n)  # 各细胞周期进度 [0,1)
    divisions = [0]  # 记录分裂事件
    
    # 凋亡信号
    apop_signal = np.zeros(n)
    
    alpha, beta = 2.0, 1.0
    gamma, dec = 0.2, 0.15
    K_base = 5.0
    
    for step in range(max_steps):
        t = step * dt
        
        # 时钟推进
        cell_cycle += dt * p["timer_speed"]
        for i in range(n):
            if cell_cycle[i] >= 1.0:
                cell_cycle[i] = 0.0
                divisions.append(step)
        
        # L6: 染色质稳定性因子 (影响Notch信号)
        chrom_factor = min(1.5, max(0.5, p["memory_stab"] * 0.3 + p["heritability"] * 0.5))
        # L4: 时钟精度 (检查点是否严格)
        tim_accuracy = min(1.5, max(0.5, p["timer_precis"] * 0.5 + p["check_g1"] * 0.3))
        
        # Notch信号
        for i in range(n):
            # 邻居平均配体
            neighbors = [j for j in range(n) if j != i]
            if neighbors:
                avg_G = np.mean([G[j] * p.get("lat_inhib", 1.0) for j in neighbors])
            else:
                avg_G = 0
            raw_sig = avg_G * p["notch_sens"]
            # L6影响信号放大速率
            activation = alpha * raw_sig * (1 - N[i]) * chrom_factor
            # L5影响噪声幅度
            noise = 0.001 * (rng.rand() - 0.5) * (1.0 / max(0.1, p.get("ach_sens", 0.5))) * tim_accuracy
            N[i] += dt * (p["lin12_expr"] * activation - p.get("notch_decay", 0.5) * N[i] + noise)
            N[i] = max(0, min(2, N[i]))
        
        # L5: 神经递质对配体能力的调节
        driver_factor = min(1.3, max(0.3, p["ach_sens"] * 0.3 + p["gaba_sens"] * 0.3 + p["glu_sens"] * 0.2))
        
        # L5: 神经递质对配体能力的调节
        driver_factor = min(1.3, max(0.3, p["ach_sens"] * 0.3 + p["gaba_sens"] * 0.3 + p["glu_sens"] * 0.2))
        
        # 配体更新
        for i in range(n):
            inhibition = 1.0 + K_base * N[i]**2 / (N[i]**2 + 0.1)
            G[i] += dt * (p["lag2_expr"] * beta / inhibition - gamma * G[i] * driver_factor)
            G[i] = max(0, min(5, G[i]))
        
        # L7: 凋亡阈值
        apop_bias = min(2.0, max(0.5, p["cell_survive"] * 0.5 + p["caspase"] * 0.3))
        survival_thresh = 0.05 * apop_bias
        
        # 凋亡信号 (n_cells=3时: 半数致死)
        if n == 3:
            target_death = int(n / 2)
            # 结合Notch水平和细胞存活率的死亡判定
            death_scores = -N * p.get("cell_survive", 0.8) + G * 0.1 + p.get("caspase", 0.5) * 0.05
            death_candidates = np.argsort(-death_scores)
            for i, death_idx in enumerate(death_candidates):
                if i < target_death and N[death_idx] > -0.5 and G[death_idx] < survival_thresh * 10:
                    N[death_idx] = -1  # 标记为死亡
                    G[death_idx] = 0
            # 只有活跃细胞继续
            alive = N >= 0
            if sum(alive) <= 1 and step > 20:
                break
        
        # 收敛检测
        if step > 40:
            dN = np.max(np.abs(np.diff([np.mean(N[:n//2]), np.mean(N[n//2:])])))
            if dN < 0.001 and np.std(N) > 0.1:
                break
    
    # 结果判定
    active = N >= 0
    n_alive = sum(active)
    
    if n_alive >= 2:
        max_idx = np.argmax(N)
        min_idx = np.argmin(N) if n_alive >= 2 else max_idx
        correct = N[max_idx] > N[min_idx] + 0.05
    else:
        correct = False
    
    return {
        "correct": correct,
        "n_alive": n_alive,
        "N_final": list(round(v, 3) for v in N),
        "steps": step + 1,
        "divisions": len(divisions),
        "stdev": float(np.std([v for v in N if v >= 0])) if sum(active) > 0 else 0,
    }


# ── 多任务适应性 ────────────────────────────────────────

def fitness_v04(network, n_seeds=5, verbose=False):
    """多任务适应性函数。"""
    
    total = 0.0
    log = {}
    
    # 任务1: AC/VU 双细胞 (权重1.0)
    scores_2cell = []
    for seed in range(n_seeds):
        r = simulate_multi_cell(network, n_cells=2, seed=seed)
        if r["correct"] and r["n_alive"] >= 2:
            s = 1.0 + min(0.3, r["stdev"] * 2)
            if r["steps"] < 100:
                s += 0.1
            elif r["steps"] > 300:
                s -= 0.05
        else:
            s = 0.15 if r["n_alive"] >= 2 else 0.0
        scores_2cell.append(s)
    task1 = sum(scores_2cell) / len(scores_2cell)
    total += task1 * 1.0
    log["ac_vu"] = round(task1, 3)
    
    # 任务2: 三细胞竞争 (权重0.8)
    scores_3cell = []
    for seed in range(n_seeds):
        r = simulate_multi_cell(network, n_cells=3, seed=seed + 100)
        if r["correct"] and r["n_alive"] >= 1:
            s = 0.5 + min(0.3, r["stdev"] * 2)
        elif r["n_alive"] >= 1 and not r["correct"]:
            s = 0.15
        else:
            s = 0.0
        scores_3cell.append(s)
    task2 = sum(scores_3cell) / len(scores_3cell)
    total += task2 * 0.8
    log["three_cell"] = round(task2, 3)
    
    # 任务3: 发育速度 (0-0.3)
    r1 = simulate_multi_cell(network, n_cells=2, seed=0)
    speed = 0.3 if r1["steps"] < 80 else (0.15 if r1["steps"] < 200 else (0.0 if not r1["correct"] else 0.05))
    total += speed
    log["speed"] = round(speed, 3)
    
    # 任务4: 噪声容忍 (多种子AC/VU方差 0-0.2)
    noise_scores = []
    for seed in range(10):
        r = simulate_multi_cell(network, n_cells=2, seed=seed * 7 + 3)
        noise_scores.append(1.0 if r["correct"] and r["n_alive"] >= 2 else 0.0)
    noise_rate = sum(noise_scores) / len(noise_scores) if noise_scores else 0
    noise_score = noise_rate * 0.2
    total += noise_score
    log["noise"] = round(noise_score, 3)
    
    # 能量成本 (扣分)
    energy_cost = 0.0
    for k, v in network.params.items():
        base = WILD_PARAMS[k]
        cost = abs(v - base) * 0.02
        energy_cost += cost
    total -= energy_cost
    log["energy"] = round(-energy_cost, 3)
    
    if verbose:
        print(f"  tasks: {log}, total={total:.3f}")
    
    return max(0.0, total)


# ── 种群 + 轨迹追踪 ────────────────────────────────────

class PopulationSnapshot:
    def __init__(self, generation, individuals, scores):
        self.gen = generation
        self.n = len(individuals)
        best_idx = np.argmax(scores)
        self.best_score = scores[best_idx]
        self.avg_score = float(np.mean(scores))
        self.median = float(np.median(scores))
        self.worst = float(np.min(scores))
        self.std = float(np.std(scores))
        
        # 参数分布 (按组)
        self.param_means = {}
        self.best_params_diff = {}
        
        if individuals:
            best_ind = individuals[best_idx]
            self.best_params_diff = {
                k: round(best_ind.params[k] - WILD_PARAMS[k], 3)
                for k in WILD_PARAMS
                if abs(best_ind.params[k] - WILD_PARAMS[k]) > 0.1
            }
            
            for k in WILD_PARAMS:
                vals = [ind.params[k] for ind in individuals]
                self.param_means[k] = round(float(np.mean(vals)), 3)
    
    def to_dict(self):
        return {
            "gen": self.gen,
            "best": round(self.best_score, 3),
            "avg": round(self.avg_score, 3),
            "median": round(self.median, 3),
            "worst": round(self.worst, 3),
            "std": round(self.std, 3),
            "best_params": self.best_params_diff,
            "param_means": self.param_means,
        }


def create_population(pop_size=50):
    pop = []
    for _ in range(pop_size):
        ind = GeneNetwork()
        for k in ind.params:
            ind.params[k] *= (1 + np.random.normal(0, 0.1))
            ind.params[k] = max(0.01, min(10.0, ind.params[k]))
        pop.append(ind)
    return pop


# ── 运行一次演化 ────────────────────────────────────────

def run_evolution(pop_size=50, n_generations=200, mutation_rate=0.12,
                  severity=0.15, elitism=0.1, n_seeds=5,
                  label="test", checkpoint_every=50):
    
    pop = create_population(pop_size)
    
    history = []
    trajectory = []
    
    for gen in range(n_generations):
        scores = [fitness_v04(ind, n_seeds=n_seeds) for ind in pop]
        
        snap = PopulationSnapshot(gen, pop, scores)
        history.append(snap)
        
        if gen % 20 == 0 or gen == n_generations - 1:
            best_params = "; ".join(f"{k}={v}" for k, v in snap.best_params_diff.items())[:60]
            print(f"  Gen {gen:3d} | b={snap.best_score:.3f} a={snap.avg_score:.3f} "
                  f"w={snap.worst:.3f} | {best_params}")
        
        if gen % 10 == 0:
            trajectory.append(snap.to_dict())
        
        # 选择
        ranked_idx = sorted(range(len(scores)), key=lambda i: -scores[i])
        ranked = [(scores[i], pop[i]) for i in ranked_idx]
        
        n_elite = max(1, int(pop_size * elitism))
        new_pop = [ranked[i][1] for i in range(n_elite)]
        while len(new_pop) < pop_size:
            pool = random.choices(ranked, k=3)
            p1 = max(pool[:2], key=lambda x: x[0])[1]
            p2 = max(pool[1:3], key=lambda x: x[0])[1]
            child = p1.crossover(p2)
            child = child.mutate(mutation_rate, severity)
            new_pop.append(child)
        pop = new_pop
        
        if (gen + 1) % checkpoint_every == 0:
            data = {
                "label": label,
                "gen": gen,
                "trajectory": trajectory,
                "final_pop": [
                    {k: round(v, 3) for k, v in ind.params.items()}
                    for ind in pop[:10]
                ]
            }
            path = RESULTS_DIR / f"checkpoint_{label}_{gen:04d}.json"
            with open(path, 'w') as f:
                json.dump(data, f, indent=2)
    
    return {
        "history": history,
        "trajectory": trajectory,
        "final_pop": pop,
    }


# ── 主运行 ──────────────────────────────────────────────

if __name__ == "__main__":
    print(f"{'='*65}")
    print(f"  C. elegans OS v0.4 — 30参数 + 多任务适应性")
    print(f"  参数数: {len(WILD_PARAMS)}")
    print(f"  任务数: AC/VU(1.0) + 3细胞(0.8) + 速度(0.3) + 噪声(0.2) - 能量")
    print(f"{'='*65}")
    
    # 先跑50代验证框架
    print(f"\n  🧪 验证运行: 50代, 30参数")
    start = time.time()
    
    result = run_evolution(
        pop_size=30,
        n_generations=50,
        mutation_rate=0.12,
        n_seeds=3,
        label="v04_test"
    )
    
    elapsed = time.time() - start
    print(f"\n  ✅ 完成: {elapsed:.1f}s")
    print(f"  best={result['history'][-1].best_score:.3f} "
          f"avg={result['history'][-1].avg_score:.3f}")
    
    # 保存
    final_params = result["final_pop"][0].params if result["final_pop"] else {}
    summary = {
        "version": "0.4",
        "generations": 50,
        "pop_size": 30,
        "time_s": round(elapsed, 1),
        "final_best": result["history"][-1].best_score,
        "final_avg": result["history"][-1].avg_score,
        "best_params": {k: round(v, 3) for k, v in final_params.items()
                       if abs(v - WILD_PARAMS[k]) > 0.1},
    }
    
    path = RESULTS_DIR / "v04_test_summary.json"
    with open(path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  保存: {path}")

    # 如果验证通过, 跑全量200代
    print(f"\n  =============================================")
    print(f"  验证通过, 开始正式200代运行")
    print(f"  ============================================\n")
    
    result2 = run_evolution(
        pop_size=50,
        n_generations=200,
        mutation_rate=0.12,
        n_seeds=5,
        label="v04_full"
    )
    
    elapsed2 = time.time() - start
    print(f"\n  ✅ 全部完成: {elapsed2:.1f}s")
    print(f"  best={result2['history'][-1].best_score:.3f} "
          f"avg={result2['history'][-1].avg_score:.3f}")
    
    summary2 = {
        "version": "0.4",
        "generations": 200,
        "pop_size": 50,
        "time_s": round(elapsed2, 1),
        "final_best": result2["history"][-1].best_score,
        "final_avg": result2["history"][-1].avg_score,
        "final_worst": result2["history"][-1].worst,
        "final_std": result2["history"][-1].std,
        "best_params": {k: round(v, 3) for k, v in result2["final_pop"][0].params.items()
                       if abs(v - WILD_PARAMS[k]) > 0.1},
    }
    
    path2 = RESULTS_DIR / "v04_full_summary.json"
    with open(path2, 'w') as f:
        json.dump(summary2, f, indent=2)
    print(f"  保存: {path2}")
