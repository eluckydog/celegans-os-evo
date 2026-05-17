#!/usr/bin/env python3
"""
LLM-as-超算适配器 v0.1
把演化模拟器的轨迹封装成prompt → 喂给LLM → 解析推荐 → 输出配置

用法:
  python llm_advisor.py --report results_v04/v04c_baseline_report.json

输出:
  - LLM建议（stdout）
  - 下一步配置（JSON）
"""
import json, sys, os
from pathlib import Path

# 子任务评分 (基线 - v0.4c baseline 99代)
# 从报告文件读, 如果没有的话用硬编码的默认值

COUPLING_RULES_STR = """1. lag2_expr ↔ lin12_expr (正向, 强度0.05) — 配体-受体共演化
2. notch_decay ↔ notch_sens (正向, 0.03) — 高敏感需慢衰减
3. fork_rate ↔ timer_speed (正向, 0.04) — 分裂与时钟同步
4. cell_survive ↔ apop_thresh (正向, 0.03) — 存活高则凋亡阈值高
5. caspase ↔ egl1_sens (正向, 0.03) — 凋亡正反馈
6. h3k4me3 ↔ h3k27me3 (反向, -0.03) — 激活-抑制反相关
7. memory_stab ↔ heritability (正向, 0.03) — 记忆与可遗传
8. check_g1 ↔ check_g2 (正向, 0.03) — 检查点协同
9. mother_size ↔ daughter_asym (正向, 0.02) — 母大则不对称可大
10. ach_sens ↔ gaba_sens (正向, 0.02) — 神经递质平衡"""


def load_report(path):
    with open(path, 'r') as f:
        return json.load(f)


def build_prompt(report, title="演化轨迹分析"):
    """从报告构建prompt。"""
    prompt = f"""你正在分析一个30参数的C. elegans操作系统演化模型。

==============================
{title}
==============================

当前状态:
- 方向: {report.get("label", "baseline")}
- 种群: {report.get("gens", "?")}代
- 最终best={report['best']:.3f}, avg={report['avg']:.3f}
- 转折点: {report.get("tips", 0)}个

最佳个体参数偏移 (vs 野生型):
"""
    for k, v in report.get("best_diff", {}).items():
        prompt += f"  {k}: {'+' if v >= 0 else ''}{v:.3f}\n"

    prompt += f"""
跨层耦合规则:
{COUPLING_RULES_STR}

注意:
- 正偏移 = 参数值高于野生型 (更活跃)
- 负偏移 = 参数值低于野生型 (更抑制)
- 当前最佳个体适应度 = 子任务加权和 - 跨层耦合惩罚

请你(作为C. elegans发育生物学专家 + 演化博弈理论家):

1. **诊断**: 当前策略有什么问题? (比如“过度依赖凋亡补偿”)
2. **推荐1-2个调整方向**: 具体说调哪个参数、朝哪个方向调、为什么
3. **反直觉策略**: 基于线虫真实生物学, 有没有一个出人意料的策略值得尝试?
4. **约束建议**: 哪条跨层耦合规则可能设错了?

只输出2-3个简短推荐。"""
    return prompt


def build_strategic_prompt(history, label="next_direction"):
    """从多次迭代历史构建战略prompt。"""
    prompt = f"""你是一个赛博演化生物学家。下面是C. elegans演化模拟器最近几次迭代的表现。

"""
    for i, h in enumerate(history[-3:]):  # 最近3次
        prompt += f"\n--- 迭代{i+1}: {h['batch']}代, {h['direction']} ---\n"
        prompt += f"  起止: {h['start_best']:.2f} → {h['end_best']:.2f}\n"
        prompt += f"  关键变化: {h.get('key_changes', '无')}\n"
        prompt += f"  问题: {h.get('problem', '无')}\n"

    prompt += """
基于以上历史, 回答:
1. 这个系统在收敛还是发散?
2. 下一步应该:
   a) 保持当前方向再跑200代
   b) 换方向探索(换什么?)
   c) 加约束/改权重
3. 有没有出现你需要关注的"意外模式"?
"""
    return prompt


def parse_llm_response(response):
    """从LLM回复提取推荐。(初步实现)"""
    recommendations = []
    lines = response.strip().split('\n')
    for line in lines:
        line = line.strip()
        if line.startswith(('- ', '1.', '2.', '3.')):
            recommendations.append(line)
    if not recommendations:
        recommendations.append(response[:200])
    return recommendations


def print_prompt(prompt):
    """输出prompt到stdout。"""
    print("\n" + "="*65)
    print("  LLM PROMPT (发给当前模型)")
    print("="*65)
    print(prompt)
    print("="*65)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--report', type=str, help='Path to report JSON')
    parser.add_argument('--print', action='store_true', help='Just print prompt')
    args = parser.parse_args()

    if args.report:
        report = load_report(args.report)
        prompt = build_prompt(report)
        print_prompt(prompt)
        print("\n[要把它发给LLM, 运行 python llm_advisor.py --report <path> 得到prompt后, 手动copy到对话中]")
        print("[或者直接在这里回复LLM建议]")
    else:
        # 示例 - baseline
        rpt_path = Path(__file__).parent / "results_v04" / "v04c_baseline_report.json"
        if rpt_path.exists():
            report = load_report(rpt_path)
            prompt = build_prompt(report)
            print_prompt(prompt)
