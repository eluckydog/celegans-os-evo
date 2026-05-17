# celegans-os-evo

> **Treating *C. elegans* as a computer that has been running the same code for 450 million years — reverse-engineering its developmental operating system through evolutionary simulation.**

---

## What is this?

*C. elegans* (the roundworm) is the most thoroughly mapped multicellular organism in biology: from 1 fertilized egg to 959 cells, every division trace, apoptosis timing, and neural connection is known.

This project approaches it from two complementary angles:

### Angle 1: Reverse Software Engineering (OS Decompilation)

Map genome → development → nervous system → behavior into a layered computer operating system:

| OS Layer | Biological Analog |
|----------|------------------|
| **Bootloader** | Maternal factors (SKN-1, PIE-1) |
| **Process Tree** | Fork/apoptosis lineage from 1 cell → 959 cells |
| **IPC Stack** | Notch/Wnt/EGF/TGF-β signaling pathways |
| **System Calls** | 10 core developmental instructions |
| **Scheduler** | Cell division timing + dauer hibernation interrupt |
| **Device Drivers** | 8 neurotransmitter protocols |
| **Memory Management** | Chromatin modification + piRNA |
| **Bug Tracker** | Known mutants × OS defect reports |

> **Key findings**: AC/VU cell fate = distributed election algorithm (Bully), 131 apoptosis events = OS-level `kill()` calls, dauer larva = system S3 sleep state.

### Angle 2: Evolutionary Simulator

Treat key signaling pathway parameters as "evolvable variables" and let populations undergo natural selection in a virtual environment, observing which parameter combinations reproduce wild-type-like developmental behavior.

**30-parameter model** across 6 subsystems:
- **AC/VU** — Notch Delta lateral inhibition (two-cell election)
- **Fork** — Cell division rate & precision
- **Timer** — Developmental clock & checkpoints
- **Noise** — Neurotransmitter robustness
- **Memory** — Chromatin epigenetic memory
- **Apoptosis** — Threshold & execution

---

## Quick Start

```bash
# Requirements: Python 3.9+, numpy
git clone https://github.com/YOUR_USER/celegans-os-evo
cd celegans-os-evo

# Run baseline evolution simulation (30 params, 50 individuals, 100 generations)
python evo_sim_v04c_final.py
```

Outputs to `results_v04/`: fitness trajectory per 10 gens, tipping point detection, best-individual parameter offsets.

### Run your own direction

Edit the `directions` dict at the bottom of `evo_sim_v04c_final.py`:

```python
directions = {
    "my_exp": {"n_gens": 500, "pop_size": 30, "mut_rate": 0.15},
}
```

---

## Results Snapshot (Rounds 2-3 Best Solutions)

| Version | Constraint | Best Fitness | Key Behavior |
|---------|-----------|-------------|--------------|
| baseline (unconstrained) | none | 5.18 | egl1_sens/cell_survive runaway to 10× wild-type |
| D1 soft apoptosis cap | egl1>3 & cell>3 → -1.0 | 2.95 | cap works, but compensates to timer_speed |
| **D1v2 tighter cap** | threshold 2.0 + global max 8.0 | **4.25** | cell_survive drops to 4.3×, more balanced distribution |
| D4 balance penalty | high parameter std → penalty | 4.00 | good convergence but cell_survive still 6.8× |

**D1v2 is current best config**: apoptosis constraints + global parameter cap produce the most biologically plausible parameter distributions.

---

## Experiment Roadmap

### Completed runs

1. **D1 double knockout** (9-param): system compensates ✓ — lin-12 auto-recovers by gen 40
2. **D2 co-evolution**: two-population competition — no extra pressure under current fitness function
3. **D3 environmental oscillation**: high→low noise switching — stable fitness, no environmental memory
4. **D4 low mutation rate**: best convergence quality — validates real C. elegans low-mutation selection pressure
5. **D5 three-cell competition**: asymmetry achieved within 2 gens — two-cell is a straightforward extension
6. **Mutation rate sweep**: 0.05 / 0.12 / 0.20 — low mutation best
7. **Population size**: 20 vs 50 — small pop has more tipping points (16 vs 9)

### Whitebox Tracking

The layer-1 whitebox analyzer (`layer1_visualize.py`) produces:
- Fitness heatmaps (direction × generation)
- Tipping point detection (adaptive jumps)
- Parameter trajectory CSVs

---

## Project Structure

```
celegans-os-evo/
├── evo_sim_v04c_final.py          # Main simulator (30 params, 6 tasks, 10 coupling rules)
├── llm_advisor.py                  # LLM analysis adapter (optional)
├── layer1_visualize.py             # Whitebox visualization
├── evolution_plan.md               # Experiment design docs
├── whitebox_analysis_v1.md         # Whitebox analysis report
│
├── results/                        # Round 1 results (9 params, 5 directions)
├── results_v04/                    # Rounds 2-3 results (30 params, 4+5 directions)
│
├── knowledge_base/
│   ├── README.md
│   └── bioinf/                     # 7-layer OS decompilation
│       ├── celegans_os_v1.0.md     # OS architecture v1.0
│       ├── celegans_os_architecture.md  # Detailed architecture
│       ├── CELEGANS_REVERSE_ENGINEERING.md # Reverse engineering overview
│       ├── DELIVERY_INDEX.md       # Deliverables index
│       ├── layer1_proc_tree_v2.py  # Process tree
│       ├── layer2_ipc_stack.py     # IPC stack
│       ├── layer3_syscalls.py      # System calls
│       ├── layer4_scheduler.py     # Scheduler
│       ├── layer5_from_genes.py    # Device drivers
│       ├── layer6_memory_manager.py # Memory management
│       ├── layer7_bug_tracker.py   # Bug tracker
│       ├── ac_vu_decompile.py      # AC/VU election algorithm decompile
│       └── dna_as_programming_language.py # DNA programming language framework
```

---

## Tech Stack

- **Python 3** — pure implementation, zero external API dependencies
- **Evolutionary algorithm**: tournament selection (k=3) + elitism (10%) + Gaussian/step mutation
- **Fitness model**: 6-subsystem weighted sum — cross-layer coupling penalty
- **Runtime**: 100 gen × 50 ind ≈ 30 seconds (single core)
- **Validation**: all gene data from NCBI Entrez public API

---

## License

MIT
