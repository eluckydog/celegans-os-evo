#!/usr/bin/env python3
"""
运行器: 分页迭代演化
用法:
    python run_evolution.py --new 100     # 新建种群, 跑100代
    python run_evolution.py --continue checkpoint_gen_0050.json 100  # 从续, 再100代
"""

import sys
import json
import time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import evo_sim_v02 as evo

SAVE_DIR = Path(__file__).parent / "checkpoints"
SAVE_DIR.mkdir(exist_ok=True)

def run_new(generations=100, pop_size=50):
    print(f"\n{'='*60}")
    print(f"新建种群: {pop_size}个个体, {generations}代")
    print(f"{'='*60}")
    
    pop = evo.create_initial_population(pop_size)
    start = time.time()
    
    history, pop = evo.evolve_generations(
        pop,
        n_generations=generations,
        pop_size=pop_size,
        mutation_rate=0.12,
        severity=0.15,
        elitism=0.1,
        n_seeds=3,
        start_gen=0,
        checkpoint_every=50
    )
    
    elapsed = time.time() - start
    print(f"\n  ✅ 完成 {generations} 代: {elapsed:.1f}s")
    
    # 最终总结
    print(f"\n  === 最终适应性 ===")
    for h in history:
        if h["gen"] % 10 == 0 or h["gen"] == generations - 1:
            print(f"  Gen {h['gen']:3d}: best={h['best']:.3f} avg={h['avg']:.3f} {h['fates']}")
    
    return pop

def run_continue(checkpoint_file, generations=100):
    print(f"\n{'='*60}")
    print(f"从 {checkpoint_file} 续跑, 再 {generations} 代")
    print(f"{'='*60}")
    
    pop, start_gen, _ = evo.load_checkpoint(checkpoint_file)
    print(f"  加载: 第{start_gen}代, {len(pop)}个个体")
    
    start = time.time()
    
    history, pop = evo.evolve_generations(
        pop,
        n_generations=generations,
        pop_size=len(pop),
        mutation_rate=0.12,
        severity=0.15,
        elitism=0.1,
        n_seeds=3,
        start_gen=start_gen + 1,
        checkpoint_every=50
    )
    
    elapsed = time.time() - start
    print(f"\n  ✅ 完成 {generations} 代: {elapsed:.1f}s")
    
    for h in history:
        if h["gen"] % 10 == 0 or h["gen"] == (start_gen + generations):
            print(f"  Gen {h['gen']:3d}: best={h['best']:.3f} avg={h['avg']:.3f} {h['fates']}")
    
    return pop


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--new", type=int, default=0,
                        help="Scratch, run N generations")
    parser.add_argument("--continue", dest="cont", nargs=2, default=None,
                        metavar=("CHECKPOINT", "GENERATIONS"),
                        help="Continue from checkpoint")
    args = parser.parse_args()
    
    if args.new:
        run_new(generations=args.new)
    elif args.cont:
        cf, gens = args.cont
        # 自动找文件名
        cf_path = SAVE_DIR / cf if not Path(cf).exists() else Path(cf)
        if not cf_path.exists():
            print(f"❌ checkpoint not found: {cf_path}")
            sys.exit(1)
        run_continue(str(cf_path), int(gens))
    else:
        print("Usage: python run_evolution.py --new 100")
        print("Usage: python run_evolution.py --continue checkpoint_gen_0050.json 100")
