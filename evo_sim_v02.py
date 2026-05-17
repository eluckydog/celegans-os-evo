#!/usr/bin/env python3
"""
C. elegans OS 演化模拟器 v0.2
- 能量惩罚 (C2) — 配体/受体过量表达浪费能量
- 多噪声测试 (C3) — 每次跑多个种子取平均
- 状态存储/恢复 — 分页迭代
"""

import numpy as np
import random
import json
import os
import time
from copy import deepcopy
from pathlib import Path
from collections import defaultdict

# ── 常量 ──────────────────────────────────────────────────────────────
SAVE_DIR = Path(__file__).parent / "checkpoints"
SAVE_DIR.mkdir(exist_ok=True)

# ── 基因网络 ──────────────────────────────────────────────────────────

WILD_PARAMS = {
    "lin12_expr_rate": 1.0,       # lin-12 受体表达速率
    "lag2_expr_rate": 1.0,         # lag-2 配体表达速率  
    "notch_sensitivity": 1.0,      # 受体对信号的敏感度
    "notch_threshold": 0.5,        # Notch信号阈值
    "lateral_inhibition_strength": 1.0,  # VU抑制邻居配体的强度
    "fork_rate": 1.0,              # cdk-1/cyb-1 分裂速率
    "fork_accuracy": 0.99,         # 分裂精度
    "timer_speed": 1.0,            # lin-4/let-7 定时器速度
    "memory_stability": 0.9,       # 命运决定的持久性
}

class GeneNetwork:
    def __init__(self, params=None):
        self.params = {k: v for k, v in WILD_PARAMS.items()}
        if params:
            self.params.update(params)

    def mutate(self, rate=0.1, severity=0.1):
        child = deepcopy(self)
        for key in child.params:
            if random.random() < rate:
                p = random.random()
                if p < 0.8:
                    delta = np.random.normal(0, severity)
                    child.params[key] *= (1 + delta)
                elif p < 0.9:
                    child.params[key] *= 0.1  # 功能丧失
                else:
                    child.params[key] *= 3.0  # 功能增益
                child.params[key] = max(0.01, min(10.0, child.params[key]))
        return child

    def energies(self):
        """能量消耗: 偏离野生型越远越耗能。"""
        cost = 0.0
        for k in ["lin12_expr_rate", "lag2_expr_rate", "notch_sensitivity",
                   "lateral_inhibition_strength"]:
            d = abs(self.params[k] - WILD_PARAMS[k])
            cost += d * 0.08  # 每偏离0.1扣0.008分
        # 特别惩罚过度配体表达 (lin-12(lf)+lag2高是模型的捷径)
        if self.params["lag2_expr_rate"] > 2.0:
            cost += (self.params["lag2_expr_rate"] - 2.0) * 0.15
        return cost


# ── ODE模拟 ──────────────────────────────────────────────────────────

def simulate_ac_vu(network1, network2, seed=42, max_steps=200, dt=0.1):
    rng = np.random.RandomState(seed)
    N1, N2 = 0.1 + rng.rand() * 0.01, 0.1 + rng.rand() * 0.01
    G1, G2 = 0.5, 0.5
    
    alpha, beta = 2.0, 1.0
    gamma, dec = 0.2, 0.15
    K_base = 5.0
    
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
                dN1, dN2 = 1, 1
                if step > 21:
                    dN1 = abs(N1 - locals().get("_N1_p", N1))
                    dN2 = abs(N2 - locals().get("_N2_p", N2))
                if dN1 < 0.001 and dN2 < 0.001:
                    break
        _N1_p, _N2_p = N1, N2
    
    fate1 = "VU" if N1 > max(N2, 0.1) else "AC"
    fate2 = "VU" if N2 > max(N1, 0.1) else "AC"
    if N1 < 0.05 and N2 < 0.05:
        fate1, fate2 = "AC", "AC"
    
    return {
        "correct": fate1 != fate2,
        "fates": f"{fate1}/{fate2}",
        "N1": N1, "N2": N2,
        "G1": G1, "G2": G2,
        "steps": step + 1
    }


# ── 适应性 (C2 能量惩罚 + C3 多种子) ──────────────────

def fitness(network1, network2=None, n_seeds=3):
    """
    适应性: network1 跟 network2 跑AC/VU。
    
    关键变更 v0.2:
    - 两个子代都是演化中的个体 (不是固定野怪)
    - 如果 network2=None, 跟自己对称跑 (测试自我不对称能力)
    """
    if network2 is None:
        network2 = network1
    
    scores = []
    
    for seed in range(n_seeds):
        # network1 在左, network2 在右
        r = simulate_ac_vu(network1, network2, seed=seed)
        
        N1, N2 = r["N1"], r["N2"]
        if r["correct"]:
            s = 1.0 + min(0.3, abs(N1 - N2) * 0.5)
            s += 0.15 if r["steps"] < 40 else (0.05 if r["steps"] < 80 else -0.05)
            # lin12太低或太高都意味着过度专业化
            if N1 < 0.1 and N2 < 0.1:
                s *= 0.3
        elif r["fates"] == "AC/AC":
            s = 0.1
        else:  # VU/VU
            s = 0.05
        if max(N1, N2) > 1.8:
            s *= 0.8
        scores.append(s)
    
    base = sum(scores) / len(scores)
    
    # C2: 能量惩罚
    base -= network1.energies()
    
    if len(scores) > 1:
        variance = np.var(scores)
        base -= variance * 0.5
    
    return max(0.0, min(2.0, base))


# ── 状态存储 ──────────────────────────────────────────────────────

def save_checkpoint(population, generation, param_penalty_strength,
                    filename=None):
    data = {
        "generation": generation,
        "param_penalty_strength": param_penalty_strength,
        "population": [
            {"params": ind.params} for ind in population
        ]
    }
    if filename is None:
        filename = SAVE_DIR / f"checkpoint_gen_{generation:04d}.json"
    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"  ✅ Checkpoint saved: {filename} ({len(population)} individuals)")
    return filename


def load_checkpoint(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        data = json.load(f)
    population = [GeneNetwork(ind["params"]) for ind in data["population"]]
    return population, data["generation"], data.get("param_penalty_strength", 0.0)


# ── 演化引擎 ──────────────────────────────────────────────────────────

def create_initial_population(pop_size=50):
    population = []
    for _ in range(pop_size):
        ind = GeneNetwork()
        for key in ind.params:
            ind.params[key] *= (1 + np.random.normal(0, 0.1))
            ind.params[key] = max(0.01, min(10.0, ind.params[key]))
        population.append(ind)
    return population


def evaluate_population(population, n_seeds=3):
    scores = []
    n_pop = len(population)
    for i in range(n_pop):
        # 随机挑一个对手 (可能跟自己)
        j = random.randint(0, n_pop - 1)
        opponent = population[j]
        scores.append(fitness(population[i], opponent, n_seeds=n_seeds))
    return scores


def evolve_generations(population, n_generations=50,
                       pop_size=50, mutation_rate=0.12, severity=0.15,
                       elitism=0.1, n_seeds=3,
                       start_gen=0, checkpoint_every=50):
    """跑 n_generations 代。"""
    
    if len(population) != pop_size:
        raise ValueError(f"pop_size mismatch: {len(population)} vs {pop_size}")
    
    history = []
    
    for gen in range(start_gen, start_gen + n_generations):
        # 评估
        scores = evaluate_population(population, n_seeds=n_seeds)
        
        ranked_idx = sorted(range(len(scores)), key=lambda i: -scores[i])
        ranked = [(scores[i], population[i]) for i in ranked_idx]
        
        best = ranked[0][0]
        avg = np.mean(scores)
        worst = ranked[-1][0]
        
        # 找到最佳个体的fates (用一个随机种子跑)
        demo = simulate_ac_vu(ranked[0][1], GeneNetwork(), seed=0)
        fates = demo["fates"]
        
        history.append({
            "gen": gen,
            "best": round(best, 4),
            "avg": round(avg, 4),
            "worst": round(worst, 4),
            "fates": fates
        })
        
        if gen % 10 == 0 or gen == start_gen + n_generations - 1:
            best_ind = ranked[0][1]
            params_str = "; ".join(f"{k}={v:.2f}" 
                  for k, v in best_ind.params.items() 
                  if abs(v - WILD_PARAMS[k]) > 0.15)
            print(f"  Gen {gen:3d} | b={best:.3f} a={avg:.3f} w={worst:.3f} "
                  f"| {fates:>5} | {params_str[:60]}")
        
        # 选择
        n_elite = max(1, int(pop_size * elitism))
        new_pop = [ranked[i][1] for i in range(n_elite)]
        
        while len(new_pop) < pop_size:
            pool = random.choices(ranked, k=3)
            parent1 = max(pool[:2], key=lambda x: x[0])[1]
            parent2 = max(pool[1:3], key=lambda x: x[0])[1]
            
            child = deepcopy(parent1)
            for key in child.params:
                if random.random() < 0.5:
                    child.params[key] = parent2.params[key]
            child = child.mutate(mutation_rate, severity)
            new_pop.append(child)
        
        population = new_pop
        
        # checkpoint
        if (gen + 1) % checkpoint_every == 0:
            save_checkpoint(population, gen, 0.0)
    
    return history, population
