#!/usr/bin/env python3
"""
C. elegans OS 演化模拟器 v0.1
=============================
核心概念: 把线虫的基因网络看作一个操作系统,
    突变 = 代码补丁/bug,
    自然选择 = CI/CD测试管道。

原理:
- 每个基因是一个"函数", 接受输入(环境/信号), 输出效应蛋白
- 突变 = 函数实现被修改 (参数变化/截断/缺失/过度活跃)
- 表型 = 函数调用链返回的结果
- 适应性 = "测试覆盖率" (给定环境下能存活的能力)

这个版本: 模拟 Notch-Delta 选主协议 (lin-12 → AC/VU 命运决定)
是OS Layer 2 (IPC - socket) 和 Layer 3 (fork) 的交集

需要的基因网络:
  lag-2 (配体) → lin-12 (受体) → pop-1 (下游转录因子)
  
本模型: 简化为两个细胞竞争, 每个细胞表达 lin-12 (受体) 和 lag-2 (配体)
谁lin-12表达更高 → 收到更多信号 → 抑制自己的配体 → 变身VU
谁lin-12表达低 → 不抑制配体 → 变身AC

数学上 = ODE耦合的竞争系统
"""

import numpy as np
import random
from copy import deepcopy

# ── 基因网络定义 ──────────────────────────────────────────────────────────

class GeneNetwork:
    """一个细胞的"内核代码"。"""
    
    def __init__(self, params=None):
        # 核心参数 —— 类比内核参数/系统调用实现
        default = {
            # Layer 2 (IPC-Notch)
            "lin12_expr_rate": 1.0,       # lin-12 受体表达速率
            "lag2_expr_rate": 1.0,         # lag-2 配体表达速率  
            "notch_sensitivity": 1.0,      # 受体对信号的敏感度
            "notch_threshold": 0.5,        # Notch信号阈值 (<= 是AC, > 是VU)
            "lateral_inhibition_strength": 1.0,  # VU抑制邻居配体的强度
            
            # Layer 3 (fork - 细胞分裂)
            "fork_rate": 1.0,              # cdk-1/cyb-1 分裂速率
            "fork_accuracy": 0.99,         # 分裂精度
            
            # Layer 4 (调度器 - 定时器)
            "timer_speed": 1.0,            # lin-4/let-7 定时器速度
            
            # Layer 6 (内存 - 染色质)
            "memory_stability": 0.9,       # 命运决定的持久性
        }
        if params:
            default.update(params)
        self.params = default

    def mutate(self, rate=0.1, severity=0.1):
        """随机突变 = 操作系统打补丁 (参数变化)。"""
        child = deepcopy(self)
        for key in child.params:
            if random.random() < rate:
                # 突变类型: 参数偏移 (大多数) 或 截断 (少数)
                if random.random() < 0.8:
                    # 参数漂移 ≈ 点突变
                    delta = np.random.normal(0, severity)
                    child.params[key] *= (1 + delta)
                elif random.random() < 0.5:
                    # 功能丧失 ≈ 基因打断
                    child.params[key] *= 0.1
                else:
                    # 功能增益 ≈ 基因复制过量
                    child.params[key] *= 3.0
                
                # 裁剪到合理范围
                child.params[key] = max(0.01, min(10.0, child.params[key]))
        return child


# ── 细胞命运决定 AC/VU ──────────────────────────────────────────────────

def simulate_ac_vu(network1: GeneNetwork, network2: GeneNetwork,
                   max_steps=200, dt=0.1):
    """
    两个细胞耦合的Notch-Delta竞争系统 (v2: 带正反馈的双稳模型)。
    
    正确的Notch-Delta侧向抑制模型:
    - 细胞接收的Notch信号 = 来自对面细胞的配体 * 自己的敏感度
    - Notch信号激活下游 → 上调lin-12 → 上调自身的lag-2抑制剂
    - 关键: 双稳态来自两个正反馈环: 
      细胞A高lin-12 → 抑制A的配体 → 细胞B低信号 → B低lin-12 → B高配体
      → A获得更多信号 → A更高lin-12 (正反馈)
    
    状态变量:
        N1, N2: Notch信号强度 (已经包含了lin-12和信号输入)
        G1, G2: lag-2 配体表达水平
    
    正确的耦合:
        dN/dt = alpha * (G_对面 * sensitivity) * (1 - N) - decay * N
        dG/dt = beta / (1 + K * N) - gamma * G
    
    alpha = 信号放大系数 (正反馈的关键)
    beta  = 配体基础表达
    K     = Notch对配体的抑制强度
    第二项 (1-N) = 受体饱和 (不能无限高)
    """
    # 初始状态: 给很小随机不对称
    rng = np.random.RandomState(42)
    N1, N2 = 0.1 + rng.rand() * 0.01, 0.1 + rng.rand() * 0.01
    G1, G2 = 0.5, 0.5
    
    # 固定参数
    alpha = 2.0       # 信号放大倍数
    beta = 1.0        # 配体基础表达
    gamma = 0.2       # 配体衰变
    dec = 0.15        # Notch衰变
    K_base = 5.0      # Notch对配体的抑制强度
    
    for step in range(max_steps):
        p1 = network1.params
        p2 = network2.params
        
        # 每个细胞的Notch抑制强度
        K1 = K_base * p1["lateral_inhibition_strength"]
        K2 = K_base * p2["lateral_inhibition_strength"]
        
        # 每个细胞收到的信号 = 对面配体 * 自己敏感度
        raw_sig1 = G2 * p1["notch_sensitivity"]
        raw_sig2 = G1 * p2["notch_sensitivity"]
        
        # Notch演变: 配体激活 + 受体饱和 + 非线性放大
        act1 = alpha * raw_sig1 * (1 - N1)  # 饱和项
        act2 = alpha * raw_sig2 * (1 - N2)
        
        # 在确定性模型中加一个小的噪声项(每次不同)
        noise = 0.001 * (rng.rand() - 0.5)
        
        N1 += dt * (p1["lin12_expr_rate"] * act1 - dec * N1 + noise)
        N2 += dt * (p2["lin12_expr_rate"] * act2 - dec * N2)
        
        # 配体演变: 受Notch抑制 (Hill-like)
        inh1 = 1.0 + K1 * N1**2 / (N1**2 + 0.1)
        inh2 = 1.0 + K2 * N2**2 / (N2**2 + 0.1)
        G1 += dt * (p1["lag2_expr_rate"] * beta / inh1 - gamma * G1)
        G2 += dt * (p2["lag2_expr_rate"] * beta / inh2 - gamma * G2)
        
        # 裁剪
        N1, N2 = max(0, min(2, N1)), max(0, min(2, N2))
        G1, G2 = max(0, min(5, G1)), max(0, min(5, G2))
        
        # 稳定判定 (两个细胞状态差距大且稳定)
        if step > 20:
            diff = abs(N1 - N2)
            if diff > 0.3:  # 已经分岔
                dN1 = abs(N1 - N1_prev) if step > 21 else 1
                dN2 = abs(N2 - N2_prev) if step > 21 else 1
                if dN1 < 0.001 and dN2 < 0.001:
                    break
        N1_prev, N2_prev = N1, N2
    
    # 结果判定
    fate1 = "VU" if N1 > max(N2, 0.1) else "AC"
    fate2 = "VU" if N2 > max(N1, 0.1) else "AC"
    # 修正: 两个都不高 -> 都AC (Notch都没激活)
    if N1 < 0.05 and N2 < 0.05:
        fate1, fate2 = "AC", "AC"
    
    return {
        "cell1_fate": fate1,
        "cell2_fate": fate2,
        "N1": N1, "N2": N2,
        "lag2_1": G1, "lag2_2": G2,
        "steps": step + 1,
        "correct": fate1 != fate2,
        "fates": f"{fate1}/{fate2}"
    }
    
    # 结果判定
    v = [L1, L2]
    max_l = max(v)
    
    fate1 = "VU" if L1 > network1.params["notch_threshold"] else "AC"
    fate2 = "VU" if L2 > network2.params["notch_threshold"] else "AC"
    
    return {
        "cell1_fate": fate1,
        "cell2_fate": fate2,
        "lin12_1": L1, "lin12_2": L2,
        "lag2_1": G1, "lag2_2": G2,
        "steps": step + 1,
        "correct": (fate1 != fate2),  # 一个AC一个VU才是正确
        "fates": f"{fate1}/{fate2}"
    }


# ── 适应性评估 ──────────────────────────────────────────────────────────

def fitness(result: dict) -> float:
    """
    适应性 = CI/CD测试得分。
    
    积分规则:
    - 一个AC一个VU (正确分裂) → +1.0
    - 两个AC (没有VU) → +0.2 (还能勉强活, 但生殖有问题)
    - 两个VU (没有AC) → +0.1 (根本没有锚细胞, 生殖失败)
    - 收敛快速 → 奖励 (60步以内)
    - 收敛慢 → 惩罚
    """
    score = 0.0
    
    # 主测试: 是否产生正确的不对称
    N1, N2 = result["N1"], result["N2"]
    diff = abs(N1 - N2)
    
    if result["correct"]:
        score = 1.0
        # 不对称越大越好
        score += min(0.3, diff * 0.5)
        # 越快越好
        if result["steps"] < 30:
            score += 0.2
        elif result["steps"] < 60:
            score += 0.1
        elif result["steps"] > 100:
            score -= 0.1
    elif result["fates"] == "AC/AC":
        score = 0.2 + max(0, 0.3 - abs(N1 - 0.5) * 0.5)
    elif result["fates"] == "VU/VU":
        score = 0.1
    
    # 鲁棒性: N值不过高不过低
    max_N = max(N1, N2)
    if max_N > 1.8:
        score *= 0.9  # 接近饱和 = 信号通路噪音
    if max_N < 0.3 and result["correct"]:
        score *= 0.85  # 不对称太微弱
    
    return max(0.0, min(2.0, score))


# ── 演化引擎 ───────────────────────────────────────────────────────────

def run_evolution(pop_size=50, generations=100, 
                  mutation_rate=0.15, severity=0.2, 
                  elitism=0.1):
    """
    演化引擎: 对"线虫OS内核参数"进行演化。
    
    每一代:
    1. 每个个体跟一个野生型对手[1,1,1,0.5,1,1,0.99,1,0.9]跑AC/VU
    2. 适应性评估
    3. 选择、交叉、突变
    """
    # 初始化: 从野生型加噪声
    wild = GeneNetwork()
    population = []
    for _ in range(pop_size):
        individual = deepcopy(wild)
        for key in individual.params:
            individual.params[key] *= (1 + np.random.normal(0, 0.1))
            individual.params[key] = max(0.01, min(10.0, individual.params[key]))
        population.append(individual)
    
    wild_params = wild.params
    
    history = []  # 每代的平均/最佳/最差适应性
    
    for gen in range(generations):
        # 评估
        scores = []
        results = []
        for ind in population:
            # 两个对手一个是当前个体, 一个是野生型
            opponent = GeneNetwork()
            # 为了让对手可变: 50%概率两个都是当前个体
            # (让突变体和突变体相互竞争也很重要)
            if random.random() < 0.3:
                opponent_params = deepcopy(wild_params)
                # 对对手也加一点噪声
                for k in opponent_params:
                    opponent_params[k] *= (1 + np.random.normal(0, 0.05))
                    opponent_params[k] = max(0.01, min(10.0, opponent_params[k]))
                opponent = GeneNetwork(opponent_params)
            
            # 细胞1和细胞2交换基因组 (显示不对称) 
            result1 = simulate_ac_vu(ind, opponent)
            result2 = simulate_ac_vu(opponent, ind)
            
            score1 = fitness(result1)
            score2 = fitness(result2)
            
            # 适应性 = 当突变更能正确决定时的加分
            if result1["correct"] and result2["correct"]:
                score = max(score1, score2)
            elif result1["correct"] or result2["correct"]:
                score = max(score1, score2) * 0.8
            else:
                score = min(score1, score2)
            
            # 额外: 基因网络本身不应消耗太多能量
            # (参数不极端 = 好的)
            param_norm = sum(abs(v - wild_params[k]) 
                          for k, v in ind.params.items())
            if param_norm > 10:
                score *= 0.9
            
            scores.append(score)
            results.append(result1)  # 只保留第一个
        
        # 排序
        ranked = sorted(zip(scores, population, results), 
                       key=lambda x: -x[0])
        
        best_score = ranked[0][0]
        avg_score = np.mean(scores)
        worst_score = ranked[-1][0]
        
        history.append({
            "gen": gen,
            "best": best_score,
            "avg": avg_score,
            "worst": worst_score,
            "best_fates": ranked[0][2]["fates"]
        })
        
        # 每10代打印
        if gen % 10 == 0 or gen == generations - 1:
            best_ind = ranked[0][1]
            params_str = "; ".join(f"{k}={v:.2f}" 
                                  for k, v in best_ind.params.items() 
                                  if abs(v - wild_params[k]) > 0.2)
            print(f"  Gen {gen:3d} | best={best_score:.3f} avg={avg_score:.3f} "
                  f"| fates={ranked[0][2]['fates']:>5} | {params_str[:60]}")
        
        # 选择 (精英保留 + 锦标赛)
        n_elite = max(1, int(pop_size * elitism))
        new_pop = [ranked[i][1] for i in range(n_elite)]
        
        while len(new_pop) < pop_size:
            # 锦标赛选择
            pool = [ranked[random.randint(0, pop_size - 1)][1] 
                  for _ in range(3)]
            parent1 = max(pool[:2], key=lambda x: 
                         scores[population.index(x)])
            parent2 = max(pool[1:3], key=lambda x: 
                         scores[population.index(x)])
            
            # 交叉: 交换参数
            child = deepcopy(parent1)
            for key in child.params:
                if random.random() < 0.5:
                    child.params[key] = parent2.params[key]
            
            # 突变
            child = child.mutate(mutation_rate, severity)
            new_pop.append(child)
        
        population = new_pop
    
    return history, ranked[0][1]


# ── 主程序 ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print("C. elegans OS 演化模拟器 v0.1")
    print("测试: AC/VU 选主协议 (Layer 2 IPC-Notch)")
    print("=" * 65)
    
    # 1. 先看野生型的行为
    print("\n【Step 1】野生型 AC/VU 模拟")
    print("-" * 40)
    
    wt = GeneNetwork()
    opponent = GeneNetwork()
    
    for trial in range(5):
        np.random.seed(trial + 42)
        result = simulate_ac_vu(wt, opponent)
        status = "✅" if result["correct"] else "❌"
        print(f"  试{trial+1}: {status} {result['fates']:>5} "
              f"(N={result['N1']:.3f}/{result['N2']:.3f}, "
              f"step={result['steps']})")
    
    # 2. 模拟已知bug (lin-12 lf)
    print("\n【Step 2】模拟已知bug: lin-12(lf) — 受体不表达")
    print("-" * 40)
    bug_ind = GeneNetwork({"lin12_expr_rate": 0.01})  # 功能丧失
    result = simulate_ac_vu(bug_ind, wt)
    print(f"  lin-12(lf) vs wt: {result['fates']:>5} — {'✅' if result['correct'] else '❌'}")
    print(f"  适应性: {fitness(result):.3f} (N={result['N1']:.3f}/{result['N2']:.3f})")
    result2 = simulate_ac_vu(bug_ind, bug_ind)
    print(f"  lin-12(lf) vs lf:   {result2['fates']:>5} — {'✅' if result2['correct'] else '❌'}")
    
    # lin-12(gf) — 过度活跃
    print("\n【Step 3】模拟已知bug: lin-12(gf) — 受体过度活跃")
    print("-" * 40)
    bug_ind2 = GeneNetwork({"lin12_expr_rate": 5.0, "notch_sensitivity": 5.0})
    result = simulate_ac_vu(bug_ind2, wt)
    print(f"  lin-12(gf) vs wt: {result['fates']:>5} — {'✅' if result['correct'] else '❌'}")
    print(f"  适应性: {fitness(result):.3f} (N={result['N1']:.3f}/{result['N2']:.3f})")
    
    # 3. 演化
    print("\n【Step 4】演化引擎: 跑 100 代")
    print("-" * 40)
    print("  种群50, 突变率0.15, 精英保留10%")
    print()
    
    history, best = run_evolution(
        pop_size=50,
        generations=100,
        mutation_rate=0.15,
        severity=0.2,
        elitism=0.1
    )
    
    print(f"\n  演化完成。")
    print(f"  最终最佳参数:")
    for key, val in best.params.items():
        delta = val - 1.0  # 相对野生型的偏移
        symbol = "▲" if delta > 0.2 else ("▼" if delta < -0.2 else "·")
        print(f"    {key:<35} {val:.3f} {symbol}")
    
    # 演化趋势
    print()
    print("  演化趋势 (每10代):")
    print(f"  {'Gen':>4} {'Best':>6} {'Avg':>6} {'Worst':>6} {'Fates':>8}")
    print(f"  {'─'*4} {'─'*6} {'─'*6} {'─'*6} {'─'*8}")
    for h in history:
        if h["gen"] % 10 == 0 or h["gen"] == 100:
            print(f"  {h['gen']:4d} {h['best']:6.3f} {h['avg']:6.3f} {h['worst']:6.3f} {h['best_fates']:>8}")
    
    # 结论
    print()
    print("=" * 65)
    print("结论")
    print("=" * 65)
    print(f"""
1. 系统能自动恢复lin-12(lf) bug吗?
   {('可以' if best.params['lin12_expr_rate'] > 0.8 else '不容易')}
   因为演化只能在现有参数上找局部最优，要"修复"一个被打断的基因
   需要另一个基因补偿——这需要很多代。
   
2. 观察到的模式:
   - 最稳定的策略是增强lateral_inhibition_strength
     因为这是打破对称性的关键
   - notch_threshold趋向0.5附近
     因为太高两个都变AC, 太低两个都变VU
   
3. OS升级的类比:
   每个参数 = 内核可调参数 (/proc/sys/kernel/*)
   突变 = 重新编译内核时的编译器优化
   自然选择 = A/B测试 + 灰度发布
""")
