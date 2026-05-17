#!/usr/bin/env python3
"""
C. elegans OS 30参数模型 + 全轨迹白盒化

## 30参数设计 (横跨6个OS层)

### Layer 2 (IPC-Notch) — 7个
1. lin12_expr_rate     — 受体表达 (原)
2. lag2_expr_rate      — 配体表达 (原)
3. notch_sensitivity   — 信号敏感度 (原)
4. notch_threshold     — 阈值 (原)
5. lateral_inhibition  — 侧向抑制 (原)
6. notch_decay         — 新增: 信号衰变速率
7. dsl_basal           — 新增: 配体基础分泌 (apx-1)

### Layer 3 (fork系统调用) — 5个
8. fork_rate           — 分裂速率 (原)
9. fork_accuracy       — 分裂对称性 (原)
10. mother_size        — 新增: 母细胞分裂前大小阈值
11. daughter_asym      — 新增: 后代不对称分配
12. polarity_strength  — 新增: PAR蛋白极化强度

### Layer 4 (调度器-定时器) — 5个
13. timer_speed        — 全局定时器 (原)
14. timer_precision    — 新增: 定时器精度 (lin-42 circadian)
15. check_point_g1     — 新增: G1检查点阈值 (cdk-4/cycD)
16. check_point_g2     — 新增: G2检查点阈值 (cdk-1/cycB)
17. dna_damage_sense   — 新增: DNA损伤感知灵敏度

### Layer 5 (驱动-神经递质) — 4个
18. ach_sensitivity    — 新增: 乙酰胆碱敏感度 (unc-17)
19. gaba_sensitivity   — 新增: GABA敏感度 (unc-49)
20. glutamate_sense    — 新增: 谷氨酸敏感度 (glr-1)
21. serotonin_balance  — 新增: 血清素平衡 (tph-1)

### Layer 6 (染色质-内存) — 5个
22. memory_stability   — 染色质修饰强度 (原, mes-2/4/6)
23. h3k4me3_rate       — 新增: 激活标记沉积速率 (set-2)
24. h3k27me3_rate      — 新增: 抑制标记沉积速率 (mes-2)
25. chromatin_noise    — 新增: 表观遗传噪声
26. heritability       — 新增: 表观标记遗传给子细胞概率

### Layer 7 (凋亡-kill系统调用) — 4个
27. apoptosis_threshold— 新增: 凋亡触发阈值 (ced-9/Egl-1)
28. caspase_activity   — 新增: 半胱天冬酶活性 (ced-3)
29. egl1_sensitivity   — 新增: 凋亡因子灵敏度
30. cell_survival      — 新增: 细胞存活因子强度

## 设计原则

1. 每个OS层至少有1个参数可以被演化
2. 新参数不是随便加的——每个对应一个已知基因/通路
3. 参数的值域统一为 0.01-10.0, 野生型在0.5-1.5之间
4. 所有参数都影响同一个适应性测试: AC/VU + 发育速度 + 细胞存活

## 适应性测试

现在不只是"能不能产生AC/VU"：
- AC/VU正确性 (weight=1.0)
- 发育速度 (weight=0.3) — fork_rate整体影响时间
- 凋亡正确 (weight=0.2) — 额外的131个凋亡事件
- 噪声容忍 (weight=0.2) — 环境波动
- 能量成本 (weight=0.3) — 参数偏离野生型

## 全轨迹白盒化 (层1)

输出:
1. 每5代: 全种群30参数分布 (0.25 -> 0.5 percentile矩阵)
2. 每代: 最佳个体参数 + 适应性 + 突变记录
3. 每代: 参数之间的相关性矩阵 (谁跟谁一起突变了)
4. 转折点检测: 适应性跳变 + 参数跳变 + 哪个参数贡献了跳变
5. 演化总结: "第X代, 参数Z偏移了Y%, 在适应性上产生了+0.XX"

第一版: 先把追踪逻辑嵌进去, 再转到热图+git-log
"""

# 目前先写框架, v0.4会包含全部30参数
# 但实际验证从15参数开始, 保证不出bug

WILD_PARAMS_30 = {
    # ---- Layer 2: IPC-Notch ----
    "lin12_expr_rate": 1.0,
    "lag2_expr_rate": 1.0,
    "notch_sensitivity": 1.0,
    "notch_threshold": 0.5,
    "lateral_inhibition": 1.0,
    "notch_decay": 0.5,
    "dsl_basal": 0.5,
    # ---- Layer 3: fork ----
    "fork_rate": 1.0,
    "fork_accuracy": 0.99,
    "mother_size": 0.8,
    "daughter_asym": 0.5,
    "polarity_strength": 0.8,
    # ---- Layer 4: timer ----
    "timer_speed": 1.0,
    "timer_precision": 0.8,
    "check_point_g1": 0.7,
    "check_point_g2": 0.7,
    "dna_damage_sense": 0.8,
    # ---- Layer 5: drivers ----
    "ach_sensitivity": 0.6,
    "gaba_sensitivity": 0.6,
    "glutamate_sense": 0.5,
    "serotonin_balance": 0.5,
    # ---- Layer 6: chromatin ----
    "memory_stability": 0.9,
    "h3k4me3_rate": 0.7,
    "h3k27me3_rate": 0.7,
    "chromatin_noise": 0.3,
    "heritability": 0.8,
    # ---- Layer 7: apoptosis ----
    "apoptosis_threshold": 0.6,
    "caspase_activity": 0.5,
    "egl1_sensitivity": 0.5,
    "cell_survival": 0.8,
}

if __name__ == "__main__":
    print(f"30参数模型: {len(WILD_PARAMS_30)}个参数")
    print(f"  Layer 2 IPC-Notch:     7个")
    print(f"  Layer 3 fork:          5个")
    print(f"  Layer 4 scheduler:     5个")
    print(f"  Layer 5 drivers:       4个")
    print(f"  Layer 6 chromatin:     5个")
    print(f"  Layer 7 apoptosis:     4个")
    print()
    print("当前步: 设计完成. 需要:")
    print("  1. 重写ODS模拟器(30参数耦合)")
    print("  2. 重写适应性函数(多任务)")
    print("  3. 合入层1追踪 + git-log生成")
