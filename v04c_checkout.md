# C. elegans OS 演化模拟器 v0.4c — checkout存档

## 设计框架
- 30参数 + 6子任务 (Notch/Fork/Timer/Noise/Memory/Apoptosis)
- 10条跨层耦合规则 + 极端值惩罚
- 白盒追踪 (每10代轨迹 + 转折点检测 + CSV输出)
- 模块化适应性函数, ODE双细胞Notch侧向抑制

## 4个方向 x 100代跑完
| 方向 | best | avg | tips | 时间 |
|------|------|-----|------|------|
| baseline (50pop, 0.12mut) | 5.18 | 4.29 | 9 | 29s |
| low_mut (0.05) | 4.79 | 4.56 | 9 | 35s |
| high_mut (0.20) | 4.03 | 2.69 | 8 | 47s |
| small_pop (20pop) | 4.89 | 3.63 | 16 | 11s |

## 未解决的问题
- 跨层约束仍不够强 (极端值涌现: cell_survive=10, egl1_sens=10)
- 子任务间负交互不充分 (高凋亡+高存活互相补偿)
- 只测试了4个方向, 非系统扫描

## 文件
- evo_sim_v04c_final.py — 主模拟器 (11.6KB)
- results_v04/v04c_*.csv — 每个方向的轨迹数据
- results_v04/v04c_*_report.json — 报告

## 下一步
- 用LLM当"超算"设计 (探索策略空间/参数关联推理)
