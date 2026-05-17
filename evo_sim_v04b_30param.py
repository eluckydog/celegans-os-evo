#!/usr/bin/env python3
"""
C. elegans OS 演化模拟器 v0.4b — 30参数, 模块化子任务

设计原则: 每个OS层对应一个独立子任务, 不强制所有参数参与所有任务。
各自演化, 但共享同一个种群和演化引擎。

子任务分配:
  Layer 2 (Notch)  → AC/VU决定 (7参数)
  Layer 3 (fork)   → 细胞分裂速度 (5参数)
  Layer 4 (timer)  → 发育同步性 (5参数)
  Layer 5 (drivers)→ 噪声容忍 (4参数)
  Layer 6 (chrom)  → 记忆持久性 (5参数)
  Layer 7 (apop)   → 单细胞存活 (4参数)

每个子任务只评估相关参数, 其他参数取中值(不扣分)。
总适应度 = 5个子任务加权和。
"""

import random, json, csv, os, time, math
import numpy as np
from pathlib import Path
from copy import deepcopy
from collections import defaultdict

ROOT = Path(__file__).parent
RESULTS_DIR = ROOT / "results_v04"
RESULTS_DIR.mkdir(exist_ok=True)

# ── 30参数定义 (不变) ─────────────────────────────────

WILD_PARAMS = {
    "lin12_expr": 1.0, "lag2_expr": 1.0, "notch_sens": 1.0,
    "notch_thresh": 0.5, "lat_inhib": 1.0, "notch_decay": 0.5,
    "dsl_basal": 0.5,
    "fork_rate": 1.0, "fork_accuracy": 0.99,
    "mother_size": 0.8, "daughter_asym": 0.5, "polarity": 0.8,
    "timer_speed": 1.0, "timer_precis": 0.8,
    "check_g1": 0.7, "check_g2": 0.7, "dna_damage": 0.8,
    "ach_sens": 0.6, "gaba_sens": 0.6, "glu_sens": 0.5,
    "serotonin": 0.5,
    "memory_stab": 0.9, "h3k4me3": 0.7, "h3k27me3": 0.7,
    "chrom_noise": 0.3, "heritability": 0.8,
    "apop_thresh": 0.6, "caspase": 0.5, "egl1_sens": 0.5,
    "cell_survive": 0.8,
}

# 每个子任务关注的参数
TASK_PARAMS = {
    "ac_vu": ["lin12_expr", "lag2_expr", "notch_sens", "notch_thresh",
              "lat_inhib", "notch_decay", "dsl_basal"],
    "fork": ["fork_rate", "fork_accuracy", "mother_size", "daughter_asym", "polarity"],
    "timer": ["timer_speed", "timer_precis", "check_g1", "check_g2", "dna_damage"],
    "noise": ["ach_sens", "gaba_sens", "glu_sens", "serotonin"],
    "memory": ["memory_stab", "h3k4me3", "h3k27me3", "chrom_noise", "heritability"],
    "apoptosis": ["apop_thresh", "caspase", "egl1_sens", "cell_survive"],
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


# ── 子任务1: AC/VU (Notch) ────────────────────────────

def task_ac_vu(params, n_seeds=3):
    """标准双细胞Notch侧向抑制。"""
    p = params
    scores = []
    for seed in range(n_seeds):
        rng = np.random.RandomState(seed)
        N1, N2 = 0.1 + rng.rand()*0.01, 0.1 + rng.rand()*0.01
        G1, G2 = 0.5, 0.5
        alpha, beta, gamma, dec = 2.0, 1.0, 0.2, 0.15
        K_base = 5.0
        dt = 0.1
        st = 0
        
        for step in range(250):
            raw1 = G2 * p["notch_sens"]
            raw2 = G1 * p["notch_sens"]
            act1 = alpha * raw1 * (1 - N1)
            act2 = alpha * raw2 * (1 - N2)
            noise = 0.001 * (rng.rand() - 0.5)
            
            N1 += dt * (p["lin12_expr"] * act1 - p["notch_decay"] * N1 + noise)
            N2 += dt * (p["lin12_expr"] * act2 - p["notch_decay"] * N2 + noise)
            N1, N2 = max(0,min(2,N1)), max(0,min(2,N2))
            
            inh1 = 1.0 + K_base * N1**2 / (N1**2 + 0.1)
            inh2 = 1.0 + K_base * N2**2 / (N2**2 + 0.1)
            G1 += dt * (p["lag2_expr"] * beta / inh1 - gamma * G1)
            G2 += dt * (p["lag2_expr"] * beta / inh2 - gamma * G2)
            G1, G2 = max(0,min(5,G1)), max(0,min(5,G2))
            st = step
            
            if step > 30 and abs(N1-N2) > 0.3 and abs(N1 - N1) < 0.001 and abs(N2 - N2) < 0.001:
                break
        
        diff = abs(N1 - N2)
        correct = diff > 0.05 and N1 >= 0 and N2 >= 0
        s = 1.0 + min(0.3, diff * 0.5) if correct else 0.0
        if st < 60: s += 0.1
        scores.append(s)
    
    return np.mean(scores)


# ── 子任务2: 分裂速度 (fork) ────────────────────────

def task_fork(params, n_seeds=3):
    p = params
    base_rate = p["fork_rate"]
    accuracy = p["fork_accuracy"]
    mother = p["mother_size"]
    asym = p["daughter_asym"]
    pol = p["polarity"]
    
    scores = []
    for seed in range(n_seeds):
        rng = np.random.RandomState(seed)
        n_divisions = 0
        total_time = 0.0
        clock = 0.0
        
        while n_divisions < 5 and total_time < 500:
            clock += base_rate * 0.2
            total_time += 1
            if clock >= 1.0:
                clock = 0.0
                n_divisions += 1
                # 分裂不对称性惩罚
                if rng.rand() > accuracy:
                    n_divisions += 0.5 - rng.rand()  # 不对称分裂浪费了半个分裂
        
        score = min(1.0, n_divisions / 5.0)  # 5次分裂得1分
        if total_time < 100: score += 0.2
        scores.append(score)
    
    return np.mean(scores)


# ── 子任务3: 时钟定时器 ──────────────────────────────

def task_timer(params, n_seeds=3):
    p = params
    scores = []
    for seed in range(n_seeds):
        rng = np.random.RandomState(seed)
        speed = p["timer_speed"]
        precis = p["timer_precis"]
        g1 = p["check_g1"]
        g2 = p["check_g2"]
        dna = p["dna_damage"]
        
        # 模拟周期
        cycles = []
        for _ in range(10):
            tick = speed
            # 检查点延迟
            if rng.rand() < (1 - precis):
                tick *= 0.5
            if rng.rand() < (1 - g1):
                tick *= 0.7
            if rng.rand() < (1 - g2):
                tick *= 0.7
            if rng.rand() < (1 - dna):
                tick *= 0.5
            cycles.append(tick)
        
        avg_cycle = np.mean(cycles)
        stdev = np.std(cycles)
        score = min(1.0, avg_cycle * 0.5)
        if stdev < 0.2: score += 0.2  # 周期一致性奖励
        if avg_cycle > 0.5: score += 0.1
        scores.append(score)
    
    return np.mean(scores)


# ── 子任务4: 噪声容忍 (drivers) ───────────────────────

def task_noise(params, n_seeds=10):
    """用Notch模型测试噪声容忍, 但参数来自drivers层。"""
    p = params
    ach = p["ach_sens"]
    gaba = p["gaba_sens"]
    glu = p["glu_sens"]
    sero = p["serotonin"]
    
    activity_level = ach * 0.4 + gaba * 0.3 + glu * 0.2 + sero * 0.1
    noise_stability = ach * 0.3 + gaba * 0.4 + glu * 0.2 + sero * 0.1
    
    # 用不同噪声测试
    rng = np.random.RandomState(42)
    noise_levels = [0.001, 0.005, 0.01, 0.05, 0.1]
    scores = []
    
    for nlevel in noise_levels:
        N = 0.5 + rng.rand() * 0.1
        for _ in range(100):
            noise = nlevel * (rng.rand() - 0.5) / max(0.1, noise_stability)
            N += activity_level * 0.05 - N * 0.1 + noise
            N = max(0, min(1, N))
        # 稳定在0.5附近=好
        stability = 1.0 - abs(N - 0.5) * 2
        scores.append(max(0, min(1, stability)))
    
    return np.mean(scores)


# ── 子任务5: 记忆持久性 (chromatin) ──────────────────

def task_memory(params, n_seeds=5):
    """染色质修饰的记忆持久性。"""
    p = params
    mem = p["memory_stab"]
    h3k4 = p["h3k4me3"]
    h3k27 = p["h3k27me3"]
    noise = p["chrom_noise"]
    herit = p["heritability"]
    
    scores = []
    for seed in range(n_seeds):
        rng = np.random.RandomState(seed)
        
        # 模拟一次"激活事件"
        state = 0.0
        # 维持率
        decay = max(0.01, 0.5 - mem * 0.4)
        rate = h3k4 * 0.3 + mem * 0.3
        repress = h3k27 * 0.3
        
        for step in range(100):
            if step < 10:
                state += rate * 0.1  # 激活期
            state -= state * decay  # 衰减
            state += rng.randn() * noise * 0.1  # 噪声
            state -= repress * 0.02  # 抑制
            state = max(0, min(1, state))
        
        # 100步后应该还有记忆
        residual = state
        score = max(0, min(1, residual * 2 + herit * 0.3))
        scores.append(score)
    
    return np.mean(scores)


# ── 子任务6: 凋亡 ─────────────────────────────────────

def task_apoptosis(params, n_seeds=5):
    p = params
    scores = []
    for seed in range(n_seeds):
        rng = np.random.RandomState(seed)
        apop = p["apop_thresh"]
        casp = p["caspase"]
        egl1 = p["egl1_sens"]
        surv = p["cell_survive"]
        
        # 模拟细胞在应激下的存活概率
        stress_levels = [0.1, 0.3, 0.5, 0.7, 0.9]
        for stress in stress_levels:
            death_prob = stress * (1.0 - surv * 0.5)
            if death_prob > apop:
                death_prob *= (1 + casp * 0.3)  # 凋亡激活
            if egl1 > 0.5:
                death_prob *= (1 + (egl1 - 0.5) * 0.5)
            survive = 1.0 - death_prob
            scores.append(survive)
    
    return np.mean(scores)


# ── 适应性: 模块化和 ─────────────────────────────────

def fitness_30(network, n_seeds=3, verbose=False):
    """6个子任务加权和。"""
    
    tasks = {
        "ac_vu": task_ac_vu,
        "fork": task_fork,
        "timer": task_timer,
        "noise": task_noise,
        "memory": task_memory,
        "apoptosis": task_apoptosis,
    }
    
    weights = {
        "ac_vu": 1.0,
        "fork": 0.6,
        "timer": 0.5,
        "noise": 0.4,
        "memory": 0.4,
        "apoptosis": 0.4,
    }
    
    total = 0.0
    scores = {}
    
    for name, fn in tasks.items():
        s = fn(network.params, n_seeds=n_seeds)
        scores[name] = round(s, 3)
        total += s * weights[name]
    
    # 能量成本 (所有参数)
    energy = sum(abs(network.params[k] - WILD_PARAMS[k]) * 0.02 for k in WILD_PARAMS)
    total -= energy
    
    if verbose:
        print(f"  tasks: {scores}, energy={-energy:.3f}, total={total:.3f}")
    
    return max(0.0, total)


# ── 演化引擎 ──────────────────────────────────────────

class PopulationSnapshot:
    def __init__(self, gen, popul, scores):
        self.gen = gen
        self.n = len(popul)
        bi = np.argmax(scores)
        self.best = scores[bi]
        self.avg = float(np.mean(scores))
        self.worst = float(np.min(scores))
        self.std = float(np.std(scores))
        
        self.best_params_diff = {}
        if popul:
            self.best_params_diff = {
                k: round(popul[bi].params[k] - WILD_PARAMS[k], 3)
                for k in WILD_PARAMS
                if abs(popul[bi].params[k] - WILD_PARAMS[k]) > 0.1
            }


def create_pop(pop_size=50):
    pop = []
    for _ in range(pop_size):
        ind = GeneNetwork()
        for k in ind.params:
            ind.params[k] *= (1 + np.random.normal(0, 0.1))
            ind.params[k] = max(0.01, min(10.0, ind.params[k]))
        pop.append(ind)
    return pop


def evolve(pop_size=50, n_gens=100, mut_rate=0.12, severity=0.15,
           elitism=0.1, n_seeds=3, label="test"):
    
    pop = create_pop(pop_size)
    trajectory = []
    
    for gen in range(n_gens):
        scores = [fitness_30(ind, n_seeds=n_seeds) for ind in pop]
        snap = PopulationSnapshot(gen, pop, scores)
        
        if gen % 10 == 0 or gen == n_gens - 1:
            bp = "; ".join(f"{k}={v}" for k, v in 
                          sorted(snap.best_params_diff.items())[:6])
            print(f"  Gen {gen:3d} | b={snap.best:.3f} a={snap.avg:.3f} w={snap.worst:.3f} | {bp}")
        
        trajectory.append({
            "gen": gen, "best": snap.best, "avg": snap.avg,
            "worst": snap.worst, "std": snap.std,
        })
        
        ranked_idx = sorted(range(len(scores)), key=lambda i: -scores[i])
        ranked = [(scores[i], pop[i]) for i in ranked_idx]
        
        n_elite = max(1, int(pop_size * elitism))
        new_pop = [ranked[i][1] for i in range(n_elite)]
        while len(new_pop) < pop_size:
            pool = random.choices(ranked, k=3)
            p1 = max(pool[:2], key=lambda x: x[0])[1]
            p2 = max(pool[1:3], key=lambda x: x[0])[1]
            child = p1.crossover(p2)
            child = child.mutate(mut_rate, severity)
            new_pop.append(child)
        pop = new_pop
    
    return pop, trajectory


if __name__ == "__main__":
    print(f"{'='*65}")
    print(f"  C. elegans OS v0.4b — 30参数模块化任务")
    print(f"  6个子任务, 5个权重分数, 1个能量成本")
    print(f"  AC/VU(1.0) + Fork(0.6) + Timer(0.5) + Noise(0.4) + Memory(0.4) + Apop(0.4)")
    print(f"{'='*65}")
    
    # 测试野生型
    print(f"\n  🧪 野生型适应性:")
    wt = GeneNetwork()
    f = fitness_30(wt, n_seeds=3, verbose=True)
    
    # 测试优化起始
    opt = GeneNetwork()
    opt.params["lag2_expr"] = 0.5  # 降低配体 (已知有利)
    f2 = fitness_30(opt, n_seeds=3, verbose=True)
    
    print(f"\n  🚀 开始演化: 50个体 x 100代")
    start = time.time()
    pop, traj = evolve(pop_size=50, n_gens=100, n_seeds=3, label="v04b")
    elapsed = time.time() - start
    
    print(f"\n  ✅ 完成: {elapsed:.1f}s")
    print(f"  最终: b={traj[-1]['best']:.3f} a={traj[-1]['avg']:.3f}")
    
    # 保存
    result = {
        "version": "0.4b",
        "pop_size": 50, "generations": 100,
        "time_s": round(elapsed, 1),
        "final_best": traj[-1]["best"],
        "final_avg": traj[-1]["avg"],
        "best_params": {k: round(pop[0].params[k] - WILD_PARAMS[k], 3)
                       for k in WILD_PARAMS
                       if abs(pop[0].params[k] - WILD_PARAMS[k]) > 0.1},
    }
    
    path = RESULTS_DIR / "v04b_summary.json"
    with open(path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  保存: {path}")
