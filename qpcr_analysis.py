# -*- coding: utf-8 -*-
"""
qPCR 数据中心 · 分析层

把一个实验（原始孔 + 样本分组 + 参数）算成可直接给界面和 Excel 用的结果字典。
计算全部复用 qpcr_core（QC、ΔΔCt、t 检验、ANOVA + Tukey、字母法）。
"""

import hashlib
import json
import math
import re
from collections import OrderedDict

import qpcr_core as q

DEFAULT_PARAMS = {
    'thr': 0.5,          # 技术重复离散阈值
    'metric': 'sd',      # sd 或 range
    'drop': True,        # 三选二剔除离群孔
    'statval': 'dcts',   # 统计所用值：dcts 或 rqs
    'comp': 'tukey',     # ≥3 组的多重比较：tukey 全部两两，dunnett 各组 vs 对照
}

COMMON_REF = ['ef', 'ef1a', 'ef1', 'eef1a', 'gapdh', 'actin', 'bactin', 'actb',
              'rpl32', 'rp49', 'rps3', 'rpl13', 'tubulin', 'tub', '18s', 'gapd',
              'ubiquitin', 'ubi']
COMMON_CTRL = ['egfp', 'gfp', 'dsegfp', 'dsgfp', 'ck', 'wt', 'control', 'ctrl',
               'con', 'mock', 'nc', 'blank', 'h2o', 'dmso', 'cont']


def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(s))]


def guess(options, keywords):
    """按词边界猜内参或对照，避免 'ef' 误命中 'Defensin' 这类内部子串。"""
    norm = lambda s: re.sub(r'[^a-z0-9]', '', str(s).lower())
    opts = [(o, norm(o)) for o in options]
    for o, on in opts:
        if on in keywords:
            return o
    for o, on in opts:
        for kw in keywords:
            if on.startswith(kw) or on.endswith(kw):
                return o
    return None


def _num(x):
    """把 NaN / inf 变成 None，保证能 JSON 序列化。"""
    if x is None:
        return None
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def ordered_groups(group_names, ctrl, group_settings):
    """分组展示顺序：有设置按设置的顺序号，否则对照组在前，其余自然排序。"""
    order = {g['group']: g.get('order') for g in (group_settings or [])}
    def key(g):
        try:
            return (0, float(order.get(g)), natural_key(g))
        except (TypeError, ValueError):
            pass
        return (1, 0 if g == ctrl else 1, natural_key(g))
    return sorted(group_names, key=key)


# 统计算法有改动时加一，让磁盘上的旧分析缓存失效
STATS_VERSION = 2

_MEMO = OrderedDict()
_MEMO_MAX = 128


def input_key(meta, wells, sample_group, group_settings=None):
    """分析输入的指纹：参数、孔、分组、分组设置不变则结果不变。"""
    m = {k: meta.get(k) for k in ('ref', 'ctrl', 'thr', 'metric', 'drop', 'statval', 'comp')}
    w = [(x.get('sample'), x.get('target'), x.get('cq'), bool(x.get('excluded')), x.get('well')) for x in wells]
    g = sorted((x.get('group'), str(x.get('order') or '')) for x in (group_settings or []))
    raw = json.dumps([STATS_VERSION, m, w, sorted(sample_group.items()), g], ensure_ascii=False, default=str)
    return hashlib.sha1(raw.encode('utf-8')).hexdigest()


def analyze_cached(meta, wells, sample_group, group_settings=None, disk=None):
    """带缓存的 analyze。disk 为可选的持久化 dict（项目级缓存）。返回 (key, result)。"""
    key = input_key(meta, wells, sample_group, group_settings)
    if key in _MEMO:
        _MEMO.move_to_end(key)
        return key, _MEMO[key]
    if disk is not None and key in disk:
        res = disk[key]
    else:
        res = analyze(meta, wells, sample_group, group_settings)
    _MEMO[key] = res
    while len(_MEMO) > _MEMO_MAX:
        _MEMO.popitem(last=False)
    return key, res


def analyze(meta, wells, sample_group, group_settings=None):
    """meta: 实验参数 dict（ref, ctrl, thr, metric, drop, statval, comp）。
    wells: list[dict]，含 target/sample/cq/excluded。
    sample_group: dict{sample: group}。
    返回结果 dict（全部可 JSON 序列化）。"""
    p = dict(DEFAULT_PARAMS)
    p.update({k: meta[k] for k in DEFAULT_PARAMS if meta.get(k) not in (None, '')})
    thr = float(p['thr'])
    drop = p['drop'] in (True, 'True', 'true', 1, '1', '是')
    ref = meta.get('ref') or ''
    ctrl = meta.get('ctrl') or ''
    warnings = []

    # 无样本名的孔（NTC 等）不参与分析
    used_wells = [w for w in wells if not w.get('excluded') and w.get('sample')]
    manual = {}
    for w in wells:
        if w.get('excluded') and w.get('cq') is not None:
            manual.setdefault((w['sample'], w['target']), []).append(w['cq'])

    groups = q.group_replicates(used_wells)
    qc_sum = q.run_qc(groups, sd_threshold=thr, metric=p['metric'], remove_outlier=drop)
    qc = []
    for g in groups:
        if g.n_raw < 2:
            status = 'single'
        elif g.qc_pass and g.outliers:
            status = 'fixed'
        elif g.qc_pass:
            status = 'pass'
        else:
            status = 'fail'
        qc.append({
            'sample': g.sample, 'target': g.target, 'group': sample_group.get(g.sample, ''),
            'n_raw': g.n_raw, 'mean': _num(g.mean) if g.cqs_used else None,
            'sd': _num(g.sd), 'range': _num(g.cq_range),
            'used': g.cqs_used, 'outliers': g.outliers,
            'manual': manual.get((g.sample, g.target), []),
            'wells': [x for x in g.wells if x],
            'status': status, 'reason': g.qc_reason,
        })

    targets_all = sorted({w['target'] for w in wells if w.get('target')}, key=natural_key)
    result = {
        'ref': ref, 'ctrl': ctrl, 'params': p,
        'targets_all': targets_all, 'targets': [], 'groups_order': [],
        'qc': qc, 'qc_summary': qc_sum, 'per_sample': [], 'genes': {},
        'warnings': warnings, 'methods_text': '',
    }

    n_empty = sum(1 for w in wells if w.get('cq') is None)
    if n_empty:
        warnings.append(f'{n_empty} 个孔无 Cq 值（无扩增或 NTC），已自动忽略。')
    ntc_amp = [w for w in wells if not w.get('sample') and w.get('cq') is not None]
    if ntc_amp:
        desc = '，'.join(f'{w.get("target")} {w.get("well")} Cq {w["cq"]:.1f}' for w in ntc_amp[:4])
        warnings.append(f'{len(ntc_amp)} 个 NTC / 无样本名孔出现扩增（{desc}），请留意污染或引物二聚体。')
    if not ref or ref not in targets_all:
        warnings.append('尚未选择有效的内参基因，无法计算 ΔΔCt。')
        return result
    if not [t for t in targets_all if t != ref]:
        warnings.append(f'内参是「{ref}」，但没有其它基因可作目的基因。是否漏导了目的基因的板？')
        return result

    dd = q.compute_ddct(groups, ref_target=ref, control_group=ctrl,
                        group_of=sample_group, sd_threshold=thr)
    n_fail = qc_sum['fail']
    if n_fail:
        warnings.append(f'{n_fail} 组技术重复离散超阈（红色），结果已纳入，建议在 QC 页检查或手动剔除孔。')
    if dd['skipped_no_ref']:
        ss = dd['skipped_no_ref']
        warnings.append(f'{len(ss)} 个样本缺内参 Cq 未纳入：{", ".join(ss[:6])}{" 等" if len(ss) > 6 else ""}')
    for t, ss in dd['skipped_no_target'].items():
        warnings.append(f'基因「{t}」有 {len(ss)} 个样本缺 Cq，未纳入该基因。')
    if dd['control_missing']:
        warnings.append(f'归一组「{ctrl}」在这些基因上没有样本，无法归一：{", ".join(dd["control_missing"])}')

    for v in dd['per_sample'].values():
        result['per_sample'].append({k: _num(v[k]) if isinstance(v[k], float) else v[k]
                                     for k in ('sample', 'group', 'target', 'ref_cq',
                                               'tgt_cq', 'dct', 'ddct', 'rq')})

    all_groups = set()
    for target in dd['targets']:
        pg = dd['per_group'][target]
        if not pg:
            continue
        all_groups.update(pg.keys())
        st = q.run_stats(pg, value_key=p['statval'], control_group=ctrl, comp=p['comp'])
        p_vs = {}
        for g in pg:
            pc = next((c for c in st['pairwise'] if {c['group_i'], c['group_j']} == {g, ctrl}), None)
            p_vs[g] = _num(pc['p']) if pc else None
        gd = {}
        for g, d in pg.items():
            ddcts = [-math.log2(r) for r in d['rqs'] if r and r > 0]
            gd[g] = {
                'n': d['n'], 'mean': _num(d['mean']), 'sem': _num(d['sem']), 'sd': _num(d['sd']),
                'rqs': d['rqs'], 'dcts': d['dcts'], 'samples': d['samples'],
                'log2fc': _num(-q.mean(ddcts)) if ddcts else None,
                'log2_points': [-x for x in ddcts],
            }
        small = [g for g, d in pg.items() if d['n'] < 2]
        result['genes'][target] = {
            'groups': gd,
            'stats': {
                'method': st['method'], 'method_label': q.METHOD_NAMES.get(st['method'], '未检验'),
                'p': _num(st['p']), 'n_groups': st['n_groups'],
                'pairwise': [{**c, 'p': _num(c['p'])} for c in st['pairwise']],
                'letters': st['letters'], 'stars_vs_control': st['stars_vs_control'],
                'p_vs_control': p_vs, 'small_groups': small,
            },
        }
        result['targets'].append(target)
    result['groups_order'] = ordered_groups(all_groups, ctrl, group_settings)
    result['methods_text'] = methods_text(p, {t['stats']['method'] for t in result['genes'].values()})
    small_any = sorted({g for t in result['genes'].values() for g in t['stats']['small_groups']})
    if small_any:
        warnings.append(f'这些分组只有 1 个生物学重复，不参与显著性检验：{", ".join(small_any)}')
    return result


def methods_text(params, methods):
    """按本实验实际用到的检验生成一段统计方法说明，供结果页显示和图注引用。"""
    val = 'RQ' if params['statval'] == 'rqs' else 'ΔCt'
    parts = []
    if 't-test' in methods:
        parts.append('两组比较采用双侧 Student t 检验')
    if 'anova' in methods:
        parts.append('多组比较采用单因素 ANOVA，组间两两比较用 Tukey HSD 校正')
    if 'dunnett' in methods:
        parts.append('多组比较采用单因素 ANOVA，各组与对照组的比较用 Dunnett 检验校正')
    if not parts:
        return ''
    return (f'以各生物学重复的 {val} 为检验值，' + '，'.join(parts) +
            '。均为参数检验，假定各组方差相等。* p<0.05，** p<0.01，*** p<0.001。')


def stats_for(values_by_group, statval_values=None, control=None):
    """对任意 {组: [值]} 做检验（时间序列逐时间点比较用）。返回 run_stats 结果。"""
    pg = {g: {'dcts': v, 'rqs': (statval_values or {}).get(g, v)} for g, v in values_by_group.items()}
    return q.run_stats(pg, value_key='dcts', control_group=control)
