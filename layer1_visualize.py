#!/usr/bin/env python3
"""
层1: 白盒追踪可视化 + Git-log生成

输入: projects/celegans_evo/results/trajectory_*.json
输出: 
  1. 演化热图 (参数分布随代际变化)
  2. 适应性轨迹图
  3. 转折点标记
  4. Git-style演化日志 (每代最关键的突变)
"""

import json, csv, os, sys
from pathlib import Path
import numpy as np

RESULTS_DIR = Path(__file__).parent / "results"

# ── 参数组 ─────────────────────────────────────────────

PARAM_GROUPS = {
    "IPC-Notch": ["lin12_expr_rate", "lag2_expr_rate", "notch_sensitivity", 
                  "notch_threshold", "lateral_inhibition_strength"],
    "fork": ["fork_rate", "fork_accuracy"],
    "timer": ["timer_speed", "memory_stability"],
}

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


def load_trajectories():
    """加载所有 trajectory_*.json 文件。"""
    trajs = {}
    for f in sorted(RESULTS_DIR.glob("trajectory_*.json")):
        with open(f, 'r') as fh:
            data = json.load(fh)
        key = data.get("direction", f.stem)
        trajs[key] = data
    return trajs


def build_git_log(trajs):
    """为每个方向生成Git-style演化日志。"""
    for key, data in trajs.items():
        label = key.replace("_", " ")
        print(f"\n{'='*65}")
        print(f"  Git Log: {label}")
        print(f"{'='*65}")
        
        trajectory = data.get("trajectory", [])
        snapshots = data.get("snapshots", [])
        
        if not trajectory:
            print("  (no trajectory data)")
            continue
        
        # 从前向后扫描每5代的snapshot
        prev_best_params = {}
        prev_avg = 0
        
        for snap in trajectory:
            gen = snap["gen"]
            best = snap["best"]
            avg = snap["avg"]
            worst = snap["worst"]
            std = snap["std"]
            params = snap.get("best_params", {})
            
            # 计算delta
            dbest = best - prev_best_params.get("best", best)
            davg = avg - prev_avg
            
            # 突变摘要
            changed_params = []
            for k, v in sorted(params.items()):
                direction = "↑" if v > 0 else "↓"
                param_name = k.replace("_", " ")
                changed_params.append(f"{param_name} {direction}{abs(v):.2f}")
            
            # 输出
            fate_line = f" | avg={avg:.3f}" if gen >= 0 else ""
            severity = ""
            if abs(dbest) >= 0.2:
                severity = " ★★★ 转折点"
            elif abs(dbest) >= 0.1:
                severity = " ★★ 显著"
            elif abs(dbest) >= 0.05:
                severity = " ★ 小幅"

            # 适应性趋势
            adapt_str = ""
            if dbest > 0.05:
                adapt_str = f" +{dbest:.3f}"
            elif dbest < -0.05:
                adapt_str = f" {dbest:.3f}"
            
            # 种群结构
            pop_str = f" [μ={avg:.3f},σ={std:.3f},w={worst:.3f}]"
            
            # 参数变化
            param_str = ""
            if changed_params:
                param_str = "  " + ", ".join(changed_params[:4])
                if len(changed_params) > 4:
                    param_str += f"... (+{len(changed_params)-4} more)"
            
            print(f"  Gen {gen:3d} | b={best:.3f}{adapt_str:>8}{pop_str}")
            if param_str:
                print(f"          {param_str}")
            
            prev_best_params = {"best": best}
            prev_avg = avg


def build_param_heatmap(trajs, output_dir):
    """为每个方向生成参数变化热图 (CSV格式, 可用Excel热图)。"""
    for key, data in trajs.items():
        trajectory = data.get("trajectory", [])
        if not trajectory:
            continue
        
        # 提取所有参数
        all_params = set()
        for snap in trajectory:
            for k in snap.get("param_dist", {}):
                all_params.add(k)
        all_params = sorted(all_params)
        
        csv_path = output_dir / f"heatmap_{key}.csv"
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            # 表头: gen, param1_mean, param1_std, param2_mean, ...
            header = ["gen"]
            for p in all_params:
                header += [f"{p}_mean", f"{p}_median", f"{p}_std"]
            w.writerow(header)
            
            for snap in trajectory:
                gen = snap["gen"]
                row = [gen]
                for p in all_params:
                    dist = snap.get("param_dist", {}).get(p, {})
                    row += [
                        round(dist.get("mean", 0), 4),
                        round(dist.get("median", 0), 4),
                        round(dist.get("std", 0), 4)
                    ]
                w.writerow(row)
        
        print(f"  热图: {csv_path.name}  ({len(all_params)}参数×{len(trajectory)}代)")
    
    return output_dir


def build_breakpoint_report(trajs, output_dir):
    """转折点报告 (每方向)。"""
    for key, data in trajs.items():
        breakpoints = data.get("breakpoints", [])
        if not breakpoints:
            continue
        
        txt_path = output_dir / f"breakpoints_{key}.txt"
        lines = [f"转折点报告: {key}", f"总数: {len(breakpoints)}", ""]
        
        for i, bp in enumerate(breakpoints):
            lines.append(f"--- 转折 {i+1} ---")
            lines.append(f"  发生: Gen {bp.get('gen_before','?')} → Gen {bp.get('gen_after','?')}")
            lines.append(f"  幅度: Δbest={bp.get('delta_best','?'):+}, Δavg={bp.get('delta_avg','?'):+}")
            
            before = bp.get("before_best_params", {})
            after = bp.get("after_best_params", {})
            
            if before and after:
                lines.append("  参数变化:")
                for k in set(list(before.keys()) + list(after.keys())):
                    diff = after.get(k, 0) - before.get(k, 0)
                    if abs(diff) > 0.01:
                        arrow = "↑" if diff > 0 else "↓"
                        lines.append(f"    {k}: {before.get(k, 0):+.3f} → {after.get(k, 0):+.3f} ({arrow}{abs(diff):.3f})")
            
            lines.append("")
        
        with open(txt_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        
        print(f"  转折报告: {txt_path.name}")


def build_adaptation_curve(trajs, output_dir):
    """适应性趋势曲线 (CSV)。"""
    csv_path = output_dir / "adaptation_curves.csv"
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["direction", "label", "gen", "best", "avg", "worst", "std"])
        
        for key, data in trajs.items():
            trajectory = data.get("trajectory", [])
            snapshots = data.get("snapshots", [])
            label = data.get("direction", key)
            
            # 用trajectory (每5代)
            for snap in trajectory:
                w.writerow([
                    key,
                    label.replace("_", " "),
                    snap["gen"],
                    snap["best"],
                    snap["avg"],
                    snap["worst"],
                    snap["std"]
                ])
    
    print(f"  适应性曲线: {csv_path.name}")


def print_final_params(trajs):
    """各方向最终最佳参数对比。"""
    print(f"\n{'='*65}")
    print(f"  各方向最终参数对比")
    print(f"{'='*65}")
    
    all_keys = sorted(WILD_PARAMS.keys())
    
    # 表头
    print(f"  {'方向':<20} {'Gen':>4} {'best':>6} ", end="")
    for k in all_keys:
        print(f"{k[:10]:>10} ", end="")
    print()
    print(f"  {'─'*20} {'─'*4} {'─'*6} ", end="")
    for _ in all_keys:
        print(f"{'─'*10} ", end="")
    print()
    
    for key, data in trajs.items():
        trajectory = data.get("trajectory", [])
        snapshots = data.get("snapshots", [])
        if not trajectory:
            continue
        
        last = trajectory[-1]
        best_params = last.get("best_params", {})
        
        label = key.replace("_", " ")
        print(f"  {label:<20} {last['gen']:>4} {last['best']:>6.3f} ", end="")
        for k in all_keys:
            v = best_params.get(k, 0)
            if abs(v) > 0.01:
                print(f"{v:>+10.3f} ", end="")
            else:
                print(f"{'·':>10} ", end="")
        print()
    
    print()
    print(f"  * 数值是相对野生型的偏移 (负=低于野生型, 正=高于野生型)")


def generate_direction_html(trajs, output_dir):
    """生成一个简短HTML报告 (简单风格, 无依赖)。"""
    html_path = output_dir / "evolution_report.html"
    
    lines = [
        '<!DOCTYPE html><html><head><meta charset="utf-8">',
        '<title>C. elegans OS 演化报告</title>',
        '<style>body{font-family:sans-serif;max-width:960px;margin:auto;padding:20px}',
        '.section{margin:20px 0;padding:15px;border:1px solid #ddd;border-radius:8px}',
        '.good{color:green}.bad{color:red}.neutral{color:gray}',
        'table{border-collapse:collapse;width:100%}',
        'td,th{border:1px solid #ddd;padding:4px 8px;text-align:right}',
        'th{background:#f5f5f5}',
        '</style></head><body>',
        '<h1>🧬 C. elegans OS 演化模拟报告</h1>',
        f'<p>生成时间: 2026-05-16</p>',
    ]
    
    # 方向摘要表
    lines.append('<div class="section"><h2>📋 方向摘要</h2><table>')
    lines.append('<tr><th>方向</th><th>代</th><th>终best</th><th>终avg</th><th>转折</th><th>结果</th></tr>')
    
    # 从phase2_summary.json加载
    summary_path = RESULTS_DIR / "phase2_summary.json"
    if summary_path.exists():
        with open(summary_path, 'r') as f:
            summaries = json.load(f)
        for s in summaries:
            bp = s.get("n_breakpoints", 0)
            result = "✅" if s["final_best"] > 1.4 else "⚠️" if s["final_best"] > 1.0 else "❌"
            lines.append(f'<tr><td style="text-align:left">{s["label"]}</td>'
                        f'<td>{s["generations"]}</td><td>{s["final_best"]:.3f}</td>'
                        f'<td>{s["final_avg"]:.3f}</td><td>{bp}</td><td>{result}</td></tr>')
    
    lines.append('</table></div>')
    
    # Git-log片段
    lines.append('<div class="section"><h2>📜 演化日志 (片段)</h2><pre style="font-size:12px">')
    for key, data in trajs.items():
        trajectory = data.get("trajectory", [])
        if not trajectory:
            continue
        label = key.replace("_", " ")
        lines.append(f'\n=== {label} ===')
        
        prev_best = 0
        for snap in trajectory[:10]:  # 只取前10个
            gen = snap["gen"]
            best = snap["best"]
            delta = best - prev_best
            arrow = "▲" if delta > 0.03 else ("▼" if delta < -0.03 else "─")
            pop_str = f"[μ={snap['avg']:.3f},σ={snap['std']:.3f}]"
            lines.append(f'Gen {gen:3d} {arrow} best={best:.3f} {pop_str}')
            prev_best = best
        
        if len(trajectory) > 10:
            lines.append('  ...')
            last = trajectory[-1]
            lines.append(f'Gen {last["gen"]:3d} best={last["best"]:.3f} [μ={last["avg"]:.3f}]')
    
    lines.append('</pre></div>')
    
    # 转折点
    lines.append('<div class="section"><h2>🔍 关键转折点</h2><table>')
    lines.append('<tr><th>方向</th><th>Gen</th><th>幅度</th><th>变化参数</th></tr>')
    for key, data in trajs.items():
        bps = data.get("breakpoints", [])
        for bp in bps[:3]:
            before = bp.get("before_best_params", {})
            after = bp.get("after_best_params", {})
            changes = []
            for k in set(list(before.keys())[:2] + list(after.keys())[:2]):
                diff = after.get(k, 0) - before.get(k, 0)
                if abs(diff) > 0.1:
                    changes.append(f"{k} {'↑' if diff>0 else '↓'}{abs(diff):.2f}")
            change_str = ", ".join(changes[:3])
            lines.append(f'<tr><td style="text-align:left">{key}</td>'
                        f'<td>{bp["gen_after"]}</td>'
                        f'<td>{"+"if bp["delta_best"]>0 else ""}{bp["delta_best"]:.3f}</td>'
                        f'<td style="text-align:left;font-size:12px">{change_str}</td></tr>')
    
    lines.append('</table></div>')
    
    lines.append('</body></html>')
    
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    
    print(f"  HTML报告: {html_path.name}")


def main():
    print(f"{'='*65}")
    print(f"  层1: 白盒追踪 + Git-log生成器")
    print(f"{'='*65}")
    
    # 加载
    trajs = load_trajectories()
    print(f"\n  加载 {len(trajs)} 个方向的轨迹\n")
    
    # 1. Git-log
    print(f"\n{'─'*65}")
    print(f"  1. Git-log 演化日志")
    print(f"{'─'*65}")
    build_git_log(trajs)
    
    # 2. 热图 (CSV)
    print(f"\n{'─'*65}")
    print(f"  2. 参数热图 (CSV)")
    print(f"{'─'*65}")
    build_param_heatmap(trajs, RESULTS_DIR)
    
    # 3. 转折点报告
    print(f"\n{'─'*65}")
    print(f"  3. 转折点报告")
    print(f"{'─'*65}")
    build_breakpoint_report(trajs, RESULTS_DIR)
    
    # 4. 适应性曲线
    print(f"\n{'─'*65}")
    print(f"  4. 适应性趋势")
    print(f"{'─'*65}")
    build_adaptation_curve(trajs, RESULTS_DIR)
    
    # 5. 最终参数对比
    print(f"\n{'─'*65}")
    print(f"  5. 最终参数对比")
    print(f"{'─'*65}")
    print_final_params(trajs)
    
    # 6. HTML报告
    print(f"\n{'─'*65}")
    print(f"  6. HTML报告")
    print(f"{'─'*65}")
    generate_direction_html(trajs, RESULTS_DIR)
    
    print(f"\n{'='*65}")
    print(f"  完成. 报告在: {RESULTS_DIR}")
    print(f"  主要产出:")
    print(f"    - evolution_report.html")
    print(f"    - adaptation_curves.csv")
    print(f"    - heatmap_*.csv (每方向)")
    print(f"    - breakpoints_*.txt (每方向)")
    print(f"{'='*65}")


if __name__ == "__main__":
    main()
