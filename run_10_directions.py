#!/usr/bin/env python3
"""
10个演化方向的批量运行器

每个方向定义:
- 哪些参数可以突变 (其他固定)
- fitness函数中的特定约束
- 预期产出 (跟哪些实际数据比对)

同时产生: 结果CSV, 突变轨迹JSON, 比对报告
"""

import sys, json, csv, os, time, random
import numpy as np
from pathlib import Path
from copy import deepcopy
sys.path.insert(0, str(Path(__file__).parent))
from evo_sim_v02 import (
    GeneNetwork, WILD_PARAMS, simulate_ac_vu,
    fitness, evaluate_population, create_initial_population,
    evolve_generations, save_checkpoint, SAVE_DIR
)

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── 10个演化方向 ────────────────────────────────────────────

# 每个方向的约束条件
DIRECTIONS = {
    "1_vu_notch_lock": {
        "label": "VU通过Notch锁定命运",
        "free_params": ["lin12_expr_rate", "lag2_expr_rate", 
                        "notch_sensitivity", "notch_threshold",
                        "lateral_inhibition_strength"],
        "fix_params": {"fork_rate": 1.0, "fork_accuracy": 0.99,
                       "timer_speed": 1.0, "memory_stability": 0.9},
        "fitness_bonus": {},  # 不加额外修饰
        "real_data": "lin-12(lf) → 2个VU; lag-2(lf) → 2个AC",
        "hypothesis": "演化应收敛到lin-12≈1.0, lag-2≈1.0, sensitivity≈1.0, 跟真实参数一致"
    },
    "2_lf_lin12_compensation": {
        "label": "lin-12功能丧失的演化补偿",
        "free_params": ["lin12_expr_rate", "lag2_expr_rate", 
                        "notch_sensitivity", "notch_threshold",
                        "lateral_inhibition_strength",
                        "fork_rate", "fork_accuracy"],  
        "fix_params": {"lin12_expr_rate": 0.01},  # 模拟lin-12打掉
        "fitness_bonus": {},  
        "real_data": "lin-12(null) 表型: 2个VU. 问题: 其他基因能否补偿?",
        "hypothesis": "lin-12被打掉 → 演化会增强lag2/notch_sensitivity来补偿"
    },
    "3_lf_lag2_compensation": {
        "label": "lag-2功能丧失的演化补偿",
        "free_params": ["lin12_expr_rate", "lag2_expr_rate",
                        "notch_sensitivity", "notch_threshold",
                        "lateral_inhibition_strength"],
        "fix_params": {"lag2_expr_rate": 0.01},  # lag-2打掉
        "fitness_bonus": {},
        "real_data": "lag-2(null) 表型: 2个AC. 补偿可能?",
        "hypothesis": "lag-2被打掉比lin-12更难补偿 (因为配体是唯一输入)"
    },
    "4_gf_lin12_escape": {
        "label": "lin-12过度活跃时的演化逃逸",
        "free_params": ["lag2_expr_rate", "notch_threshold",
                        "lateral_inhibition_strength"],
        "fix_params": {"lin12_expr_rate": 5.0, "notch_sensitivity": 5.0,
                       "fork_rate": 1.0, "fork_accuracy": 0.99,
                       "timer_speed": 1.0, "memory_stability": 0.9},
        "fitness_bonus": {"penalize_both_vu": 0.5},  # VU/VU重罚
        "real_data": "lin-12(gf) 表型: 2个VU. (lin-12 dominant)"
    },
    "5_env_noise_resistance": {
        "label": "环境噪声(温度波动)鲁棒性",
        "free_params": ["lin12_expr_rate", "lag2_expr_rate",
                        "notch_sensitivity", "notch_threshold",
                        "lateral_inhibition_strength"],
        "fix_params": {"fork_rate": 1.0, "fork_accuracy": 0.99,
                       "timer_speed": 1.0, "memory_stability": 0.9},
        "fitness_bonus": {"high_noise": True},
        "real_data": "温度敏感(ts)突变体: 25°C→异常"
    },
    "6_symmetric_self": {
        "label": "两个相同细胞能否自我不对称",
        "free_params": ["lin12_expr_rate", "lag2_expr_rate",
                        "notch_sensitivity", "notch_threshold",
                        "lateral_inhibition_strength"],
        "fix_params": {"fork_rate": 1.0, "fork_accuracy": 0.99,
                       "timer_speed": 1.0, "memory_stability": 0.9},
        "fitness_bonus": {"self_vs_self": True},
        "real_data": "WT对称→不对称. 是否纯靠噪声?"
    },
    "7_speed_vs_reliability": {
        "label": "快速决定 vs 可靠决定",
        "free_params": ["lin12_expr_rate", "lag2_expr_rate",
                        "notch_sensitivity", "notch_threshold",
                        "lateral_inhibition_strength"],
        "fix_params": {"fork_rate": 1.0, "fork_accuracy": 0.99,
                       "timer_speed": 1.0, "memory_stability": 0.9},
        "fitness_bonus": {"speed_weight": 2.0},  # 速度权重加倍
        "real_data": "实验: AC/VU决定在~2h内完成"
    },
    "8_memory_persistence": {
        "label": "命运决定后的记忆持久性",
        "free_params": ["lin12_expr_rate", "lag2_expr_rate",
                        "notch_sensitivity", "notch_threshold",
                        "lateral_inhibition_strength",
                        "memory_stability"],
        "fix_params": {"fork_rate": 1.0, "fork_accuracy": 0.99,
                       "timer_speed": 1.0},
        "fitness_bonus": {"memory_persist": True},
        "real_data": "Notch命运决定一旦锁定不反转"
    },
    "9_double_knockout": {
        "label": "双重击打: lin-12+lag-2同时丧失",
        "free_params": ["notch_sensitivity", "notch_threshold"],
        "fix_params": {"lin12_expr_rate": 0.01, "lag2_expr_rate": 0.01,
                       "fork_rate": 1.0, "fork_accuracy": 0.99,
                       "timer_speed": 1.0, "memory_stability": 0.9,
                       "lateral_inhibition_strength": 1.0},
        "fitness_bonus": {},
        "real_data": "lin-12; lag-2双突变: 无法弥补"
    },
    "10_evo_dev_link": {
        "label": "发育速率与命运判定耦合",
        "free_params": ["lin12_expr_rate", "lag2_expr_rate",
                        "notch_sensitivity", "notch_threshold",
                        "lateral_inhibition_strength",
                        "fork_rate", "fork_accuracy"],
        "fix_params": {"timer_speed": 1.0, "memory_stability": 0.9},
        "fitness_bonus": {"rate_coupling": True},
        "real_data": "温度改变→发育速率变化→决定时间窗口"
    }
}


# ── fitness包装器 (按方向微调) ──────────────────────────

def direction_fitness(network1, network2, direction, n_seeds=5):
    """基础fitness + 方向特定加分/扣分。"""
    bonus = direction["fitness_bonus"]
    base = fitness(network1, network2, n_seeds=n_seeds)
    
    if bonus.get("penalize_both_vu"):
        # 4_gf_lin12_escape: VU/VU是坏结果
        r = simulate_ac_vu(network1, network2, seed=0)
        if not r["correct"] and r["fates"] == "VU/VU":
            base *= 0.3
    
    if bonus.get("high_noise"):
        # 5_env_noise: 多个高噪声种子测试
        noisy_scores = []
        for seed in range(10, 15):
            noise_idx = seed
            r = simulate_ac_vu(network1, network2, seed=noise_idx)
            ns = 1.0 if r["correct"] else (0.1 if r["fates"] == "AC/AC" else 0.05)
            noisy_scores.append(ns)
        base = (base + sum(noisy_scores)/len(noisy_scores)) / 2
    
    if bonus.get("self_vs_self"):
        # 6_symmetric_self: 跟自己对称跑
        r = simulate_ac_vu(network1, network1, seed=0)
        if not r["correct"]:
            base *= 0.5
        # 不对称越小越好
        asym = abs(r["N1"] - r["N2"])
        base += min(0.3, asym * 0.5)
    
    if bonus.get("speed_weight"):
        # 7_speed: 速度权重加倍
        r = simulate_ac_vu(network1, network2, seed=0)
        speed_bonus = 0.3 if r["steps"] < 30 else (0.1 if r["steps"] < 60 else -0.1)
        base += speed_bonus * bonus["speed_weight"]
    
    return max(0.0, min(2.0, base))


def evaluate_direction(population, direction, n_seeds=5):
    scores = []
    for i in range(len(population)):
        j = random.randint(0, len(population) - 1)
        scores.append(direction_fitness(population[i], population[j], direction, n_seeds))
    return scores


def evolve_direction(population, direction, n_generations, pop_size=50,
                     mutation_rate=0.12, severity=0.15, elitism=0.1,
                     start_gen=0, n_seeds=5):
    history = []
    
    for gen in range(start_gen, start_gen + n_generations):
        scores = evaluate_direction(population, direction, n_seeds=n_seeds)
        
        ranked_idx = sorted(range(len(scores)), key=lambda i: -scores[i])
        ranked = [(scores[i], population[i]) for i in ranked_idx]
        
        best = ranked[0][0]
        avg = np.mean(scores)
        worst = ranked[-1][0]
        
        demo = simulate_ac_vu(ranked[0][1], GeneNetwork(), seed=0)
        
        history.append({
            "gen": gen,
            "best": round(best, 4),
            "avg": round(avg, 4),
            "worst": round(worst, 4),
            "fates": demo["fates"]
        })
        
        if gen % 10 == 0 or gen == start_gen + n_generations - 1:
            best_ind = ranked[0][1]
            changed = "; ".join(f"{k}={v:.2f}" for k, v in best_ind.params.items()
                               if abs(v - WILD_PARAMS.get(k, v)) > 0.15)
            print(f"  Gen {gen:3d} | b={best:.3f} a={avg:.3f} w={worst:.3f} | {demo['fates']:>5} | {changed[:50]}")
        
        # 选择
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
    
    return history, population


# ── 运行所有10个方向 ──────────────────────────────────

def run_all():
    summary = []
    
    for key in sorted(DIRECTIONS.keys()):
        d = DIRECTIONS[key]
        print(f"\n{'='*60}")
        print(f"方向 {key}: {d['label']}")
        print(f"  自由参数: {d['free_params']}")
        print(f"  固定参数: {d['fix_params']}")
        print(f"  bonus: {d['fitness_bonus']}")
        print(f"  实际数据: {d['real_data'][:60]}...")
        print(f"{'='*60}")
        
        # 创建初始种群, 应用fix_params
        pop = create_initial_population(50)
        for ind in pop:
            for k, v in d["fix_params"].items():
                ind.params[k] = v
        
        start = time.time()
        history, final_pop = evolve_direction(
            pop, d,
            n_generations=100,
            pop_size=50,
            n_seeds=5
        )
        elapsed = time.time() - start
        
        # 最终最佳个体
        best = simulate_ac_vu(final_pop[0], GeneNetwork(), seed=0)
        
        summary.append({
            "direction": key,
            "label": d["label"],
            "time_s": round(elapsed, 1),
            "final_best": history[-1]["best"],
            "final_avg": history[-1]["avg"],
            "final_fates": history[-1]["fates"],
            "best_params": {
                k: round(v, 3) for k, v in final_pop[0].params.items()
                if abs(v - WILD_PARAMS.get(k, v)) > 0.1
            },
            "real_data_hit": check_real_data_match(key, final_pop[0])
        })
        
        print(f"\n  结果: {elapsed:.1f}s | best={history[-1]['best']:.3f} | {history[-1]['fates']}")
        print(f"  最佳个体偏离: {', '.join(f'{k}={final_pop[0].params[k]:.2f}' for k in sorted(final_pop[0].params.keys()) if abs(final_pop[0].params[k] - WILD_PARAMS.get(k, 1.0)) > 0.15)}")
    
    # 汇总报告
    print(f"\n\n{'='*60}")
    print("10个方向汇总")
    print(f"{'='*60}")
    for s in summary:
        hit = "✅" if s["real_data_hit"] else "❌"
        print(f"  {s['direction']:>3} | {s['label']:<30} | "
              f"b={s['final_best']:.3f} | {s['final_fates']:>5} | {hit}")
    
    # 保存
    csv_path = RESULTS_DIR / "direction_summary.csv"
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(["direction", "label", "time_s", "final_best", "final_avg",
                     "final_fates", "best_params", "real_data_hit"])
        for s in summary:
            w.writerow([s["direction"], s["label"], s["time_s"],
                       s["final_best"], s["final_avg"], s["final_fates"],
                       str(s["best_params"]), s["real_data_hit"]])
    
    # JSON
    json_path = RESULTS_DIR / "direction_results.json"
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    
    print(f"\n  保存: {csv_path}")
    print(f"  保存: {json_path}")
    
    return summary


def check_real_data_match(key, best_individual):
    """检查最佳个体是否符合实际生物学数据。"""
    p = best_individual.params
    
    if key == "1_vu_notch_lock":
        # WT应该: lin12~1.0, lag2~1.0, 阈值~0.5
        if 0.5 < p["lin12_expr_rate"] < 1.5 and 0.5 < p["lag2_expr_rate"] < 1.5:
            return True
    elif key == "2_lf_lin12_compensation":
        r = simulate_ac_vu(best_individual, GeneNetwork({"lin12_expr_rate": 0.01}), seed=0)
        return r["correct"]  # 能补偿成功就算
    elif key == "3_lf_lag2_compensation":
        r = simulate_ac_vu(best_individual, GeneNetwork({"lag2_expr_rate": 0.01}), seed=0)
        return r["correct"]
    elif key == "4_gf_lin12_escape":
        r = simulate_ac_vu(best_individual, GeneNetwork({"lin12_expr_rate": 5.0, "notch_sensitivity": 5.0}), seed=0)
        return r["correct"] and "VU" not in r["fates"]  # 不能两个都VU
    elif key == "5_env_noise_resistance":
        scores = []
        for seed in range(20):
            r = simulate_ac_vu(best_individual, GeneNetwork(), seed=seed)
            scores.append(r["correct"])
        return sum(scores) >= 16  # 80%+正确率
    elif key == "6_symmetric_self":
        r = simulate_ac_vu(best_individual, best_individual, seed=0)
        return r["correct"]
    elif key == "7_speed_vs_reliability":
        r = simulate_ac_vu(best_individual, GeneNetwork(), seed=0)
        return r["steps"] < 60
    else:
        return False


if __name__ == "__main__":
    run_all()
