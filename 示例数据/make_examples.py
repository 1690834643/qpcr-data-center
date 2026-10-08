# -*- coding: utf-8 -*-
"""
生成示例数据：按 CFX Maestro「Quantification Summary」导出格式写出模拟 Cq。
数据全部为随机模拟，不对应任何真实实验。固定随机种子，重复运行结果一致。

  python 示例数据/make_examples.py
"""

import csv
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import qpcr_core as q  # noqa: E402

HEADER = ['', 'Well', 'Fluor', 'Target', 'Content', 'Sample', 'Biological Set Name',
          'Cq', 'Cq Mean', 'Cq Std. Dev', 'Starting Quantity (SQ)']
ROWS, COLS = 'ABCDEFGH', range(1, 13)


def wells():
    for r in ROWS:
        for c in COLS:
            yield f'{r}{c:02d}'


def plate(entries):
    """entries: list[(target, sample, [cq...])]，cq 为 None 表示无扩增。返回 CFX 行。"""
    it = wells()
    out = []
    for target, sample, cqs in entries:
        valid = [c for c in cqs if c is not None]
        mean = sum(valid) / len(valid) if valid else None
        sd = (sum((c - mean) ** 2 for c in valid) / (len(valid) - 1)) ** 0.5 if len(valid) > 1 else 0.0
        content = 'NTC' if sample == 'NTC' else 'Unkn'
        for c in cqs:
            out.append(['', next(it), 'SYBR', target, content, '' if sample == 'NTC' else sample, '',
                        round(c, 4) if c is not None else 'NaN',
                        round(mean, 4) if mean is not None else 'NaN', round(sd, 4), 'NaN'])
    return out


def reps(rng, mu, n=3, tech_sd=0.10):
    return [mu + rng.gauss(0, tech_sd) for _ in range(n)]


def write_xlsx(path, rows):
    q.write_xlsx(path, [('Sheet1', [HEADER] + rows)])


def write_csv(path, rows):
    with open(path, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f)
        w.writerow(HEADER)
        w.writerows(rows)


def knockdown(rng, out):
    """场景一：siRNA 敲低效率。内参板与目的基因板分开跑，演示多板合并。
    含一个离群孔（三选二自动剔除）、一组离散超阈（标红）和一个有扩增的 NTC。"""
    groups = {'siNC': {'STAT3': 0.0, 'NFKB1': 0.0},
              'siSTAT3': {'STAT3': 1.9, 'NFKB1': 0.2},
              'siNFKB1': {'STAT3': -0.3, 'NFKB1': 2.4}}
    base = {'GAPDH': 17.2, 'STAT3': 24.6, 'NFKB1': 26.1}
    ref_rows, tgt_rows = [], []
    for g, shift in groups.items():
        for i in range(1, 5):
            s = f'{g}-{i}'
            bio = rng.gauss(0, 0.25)          # 样本 cDNA 量差异，内参与目的基因同步
            ref_rows.append(('GAPDH', s, reps(rng, base['GAPDH'] + bio)))
            for t in ('STAT3', 'NFKB1'):
                cqs = reps(rng, base[t] + bio + shift[t] + rng.gauss(0, 0.3))
                if s == 'siSTAT3-2' and t == 'STAT3':
                    cqs[1] += 1.6               # 离群孔
                if s == 'siNFKB1-3' and t == 'NFKB1':
                    cqs = [cqs[0] - 0.5, cqs[1] + 0.4, cqs[2] + 1.1]   # 整组离散
                tgt_rows.append((t, s, cqs))
    ref_rows.append(('GAPDH', 'NTC', [None, None]))
    tgt_rows.append(('STAT3', 'NTC', [None, 36.8]))
    write_xlsx(os.path.join(out, '01_siRNA敲低_内参板 GAPDH -  Quantification Summary.xlsx'), plate(ref_rows))
    write_xlsx(os.path.join(out, '01_siRNA敲低_目的基因板 STAT3 NFKB1 -  Quantification Summary.xlsx'), plate(tgt_rows))


def tissue(rng, out):
    """场景二：组织表达谱，多基因，CSV 格式，分两块板。以 Kidney 为归一组。"""
    tissues = {'Kidney': {'ALB': 0.0, 'MYH7': 0.0, 'GFAP': 0.0},
               'Liver': {'ALB': -5.5, 'MYH7': 0.8, 'GFAP': 0.4},
               'Heart': {'ALB': 0.6, 'MYH7': -4.8, 'GFAP': 0.3},
               'Brain': {'ALB': 1.2, 'MYH7': 1.0, 'GFAP': -5.2},
               'Muscle': {'ALB': 0.9, 'MYH7': -2.6, 'GFAP': 0.2}}
    base = {'ACTB': 18.4, 'ALB': 27.5, 'MYH7': 28.2, 'GFAP': 27.8}
    p1, p2 = [], []
    for t, shift in tissues.items():
        for i in range(1, 4):
            s = f'{t}-{i}'
            bio = rng.gauss(0, 0.3)
            p1.append(('ACTB', s, reps(rng, base['ACTB'] + bio)))
            for gname in ('ALB', 'MYH7', 'GFAP'):
                (p1 if gname == 'ALB' else p2).append(
                    (gname, s, reps(rng, base[gname] + bio + shift[gname] + rng.gauss(0, 0.35))))
    write_csv(os.path.join(out, '02_组织表达谱_板1 ACTB ALB -  Quantification Summary.csv'), plate(p1))
    write_csv(os.path.join(out, '02_组织表达谱_板2 MYH7 GFAP -  Quantification Summary.csv'), plate(p2))


def timecourse(rng, out):
    """场景三：时间序列，Control 与 LPS 刺激 × 0/2/6/24 h，分两块板。"""
    times = [0, 2, 6, 24]
    induce = {0: 0.0, 2: -2.8, 6: -4.1, 24: -1.6}     # LPS 诱导 IL6，Cq 下降
    drift = {0: 0.0, 2: 0.1, 6: -0.2, 24: 0.1}        # 对照组随时间的小幅波动
    base = {'GAPDH': 16.9, 'IL6': 29.0}
    plates = [[], []]
    for k, t in enumerate(times):
        for trt in ('Control', 'LPS'):
            for i in range(1, 4):
                s = f'{trt}_{t}h-{i}'
                bio = rng.gauss(0, 0.25)
                p = plates[k // 2]
                p.append(('GAPDH', s, reps(rng, base['GAPDH'] + bio)))
                shift = drift[t] + (induce[t] if trt == 'LPS' else 0.0)
                p.append(('IL6', s, reps(rng, base['IL6'] + bio + shift + rng.gauss(0, 0.3))))
    write_xlsx(os.path.join(out, '03_时间序列_0-2h -  Quantification Summary.xlsx'), plate(plates[0]))
    write_xlsx(os.path.join(out, '03_时间序列_6-24h -  Quantification Summary.xlsx'), plate(plates[1]))


def main():
    rng = random.Random(20261008)
    knockdown(rng, HERE)
    tissue(rng, HERE)
    timecourse(rng, HERE)
    print('示例数据已生成：', HERE)


if __name__ == '__main__':
    main()
