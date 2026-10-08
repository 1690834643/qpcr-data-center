# -*- coding: utf-8 -*-
"""
qPCR 分析助手 —— 核心引擎（纯 Python 标准库，零第三方依赖）

包含：
  1. 统计内核   : 不完全 beta 函数 -> t 检验 / 单因素 ANOVA / Tukey HSD / 字母标记(CLD)
  2. CFX reader : 手写 xlsx 解析，兼容 CFX Maestro 导出的反斜杠/大小写 bug
  3. xlsx writer: 手写带颜色标记的 xlsx 生成
  4. 数据处理   : 技术重复分组 / QC(SD 阈值 + 三选二剔除离群) / ΔΔCt 相对定量
  5. SVG 绘图   : 期刊级柱状图(均值±SEM + 散点 + 显著性字母/星号)

本文件不 import tkinter，可在无图形界面环境下 headless 测试。
"""

import math
import re
import zipfile
import io
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape as _xml_escape

# XML 1.0 不允许的控制字符；写入前剔除，否则 Excel/openpyxl 拒绝打开
_CTRL_RE = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f]')


def _xml_text(s):
    """转义为 XML 安全文本：先剔除非法控制字符，再做实体转义。"""
    return _xml_escape(_CTRL_RE.sub('', str(s)))


# =============================================================================
# 1. 统计内核
# =============================================================================

def mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs)


def variance(xs, ddof=1):
    """样本方差，ddof=1 为无偏估计。"""
    xs = list(xs)
    n = len(xs)
    if n - ddof <= 0:
        return 0.0
    m = sum(xs) / n
    return sum((x - m) ** 2 for x in xs) / (n - ddof)


def stdev(xs, ddof=1):
    return math.sqrt(variance(xs, ddof))


def sem(xs):
    """标准误 = SD / sqrt(n)。"""
    xs = list(xs)
    n = len(xs)
    if n < 2:
        return 0.0
    return stdev(xs) / math.sqrt(n)


def _betacf(a, b, x):
    """不完全 beta 的连分数展开（Lentz 法，Numerical Recipes）。"""
    MAXIT = 300
    EPS = 3.0e-16
    FPMIN = 1.0e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < EPS:
            break
    return h


def betai(a, b, x):
    """正则化不完全 beta 函数 I_x(a, b)。"""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    bt = math.exp(lbeta + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    else:
        return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_pvalue_twosided(t, df):
    """t 分布双尾 p 值。"""
    if df <= 0:
        return float('nan')
    t = abs(t)
    x = df / (df + t * t)
    return betai(df / 2.0, 0.5, x)


def f_pvalue(F, df1, df2):
    """F 分布上尾 p 值 P(X > F)。"""
    if F <= 0:
        return 1.0
    x = df2 / (df2 + df1 * F)
    return betai(df2 / 2.0, df1 / 2.0, x)


def student_ttest(a, b):
    """Student 等方差 t 检验（合并方差）。返回 (t, df, p_twosided)。"""
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return float('nan'), float('nan'), float('nan')
    m1, m2 = mean(a), mean(b)
    v1, v2 = variance(a), variance(b)
    df = n1 + n2 - 2
    sp2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / df
    se = math.sqrt(sp2 * (1.0 / n1 + 1.0 / n2))
    if se == 0:
        if m1 == m2:
            return 0.0, float(df), 1.0
        return (float('inf') if m1 > m2 else float('-inf')), float(df), 0.0
    t = (m1 - m2) / se
    p = t_pvalue_twosided(t, df)
    return t, df, p


def oneway_anova(groups):
    """单因素方差分析。groups: list[list[float]]。
    返回 dict(F, df1, df2, p, msw, grand_mean)。"""
    groups = [list(g) for g in groups if len(g) > 0]
    k = len(groups)
    N = sum(len(g) for g in groups)
    if N == 0 or k < 2:
        return {'F': float('nan'), 'df1': max(0, k - 1), 'df2': max(0, N - k),
                'p': float('nan'), 'msw': 0.0, 'grand_mean': float('nan')}
    grand = sum(sum(g) for g in groups) / N
    ssb = sum(len(g) * (mean(g) - grand) ** 2 for g in groups)
    ssw = sum(sum((x - mean(g)) ** 2 for x in g) for g in groups)
    df1 = k - 1
    df2 = N - k
    msb = ssb / df1 if df1 > 0 else 0.0
    msw = ssw / df2 if df2 > 0 else 0.0
    if df2 <= 0:
        # 每组单点，无组内自由度，检验无意义
        return {'F': float('nan'), 'df1': df1, 'df2': df2, 'p': float('nan'),
                'msw': 0.0, 'grand_mean': grand}
    F = msb / msw if msw > 0 else float('inf')
    p = f_pvalue(F, df1, df2) if msw > 0 else 0.0
    return {'F': F, 'df1': df1, 'df2': df2, 'p': p, 'msw': msw, 'grand_mean': grand}


# ---- Tukey HSD：studentized range 分布 ----

# Gauss-Legendre 节点/权重（区间 [-1, 1]），用于数值积分
def _gauss_legendre(n):
    """返回 n 点 Gauss-Legendre 节点与权重（[-1,1]）。"""
    nodes = []
    weights = []
    m = (n + 1) // 2
    for i in range(1, m + 1):
        z = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p1 = 1.0
            p2 = 0.0
            for j in range(1, n + 1):
                p3 = p2
                p2 = p1
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j
            pp = n * (z * p1 - p2) / (z * z - 1.0)
            z1 = z
            z = z1 - p1 / pp
            if abs(z - z1) < 1e-15:
                break
        nodes.append(-z)
        nodes.append(z)
        weights.append(2.0 / ((1.0 - z * z) * pp * pp))
        weights.append(2.0 / ((1.0 - z * z) * pp * pp))
    # 排序使节点单调
    pair = sorted(zip(nodes, weights))
    return [p[0] for p in pair], [p[1] for p in pair]


_GL_NODES, _GL_WEIGHTS = _gauss_legendre(64)

_SQRT2PI = math.sqrt(2.0 * math.pi)


def _phi(x):
    """标准正态 pdf。"""
    return math.exp(-0.5 * x * x) / _SQRT2PI


def _Phi(x):
    """标准正态 cdf。"""
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def _integrate(f, lo, hi):
    """64 点 Gauss-Legendre 在 [lo, hi] 上积分。"""
    half = 0.5 * (hi - lo)
    mid = 0.5 * (hi + lo)
    s = 0.0
    for node, w in zip(_GL_NODES, _GL_WEIGHTS):
        s += w * f(mid + half * node)
    return s * half


def _prange(w, k):
    """k 个标准正态的极差 <= w 的概率 P(range <= w)。
    P = k * ∫ φ(u) [Φ(u) - Φ(u-w)]^(k-1) du 。"""
    if w <= 0:
        return 0.0

    def integrand(u):
        d = _Phi(u) - _Phi(u - w)
        if d <= 0:
            return 0.0
        return _phi(u) * d ** (k - 1)

    # 积分区间分段提高精度
    lo, hi = -8.0 - w, 8.0
    n_seg = 8
    total = 0.0
    step = (hi - lo) / n_seg
    for i in range(n_seg):
        a = lo + i * step
        b = a + step
        total += _integrate(integrand, a, b)
    return min(1.0, max(0.0, k * total))


def ptukey(q, k, df):
    """studentized range 分布 CDF：P(Q <= q)，k 组，分母自由度 df。
    对 s ~ sqrt(chi2_df / df) 积分：P = ∫ f(s) * P(range <= q*s) ds 。"""
    if q <= 0:
        return 0.0
    if q == float('inf'):
        return 1.0
    if df > 2000:  # df 很大时近似 df=∞，即极差分布本身
        return _prange(q, k)

    # s 的密度：令 c = df, s>0, f(s) = 2 (df/2)^(df/2) / Γ(df/2) * s^(df-1) exp(-df s^2 / 2)
    half_df = df / 2.0
    log_const = math.log(2.0) + half_df * math.log(half_df) - math.lgamma(half_df)

    def integrand(s):
        if s <= 0:
            return 0.0
        log_f = log_const + (df - 1) * math.log(s) - half_df * s * s
        f = math.exp(log_f)
        return f * _prange(q * s, k)

    # s 主要质量在 1 附近，积分区间 (0, ~ 1 + 几个 sd)
    hi = 1.0 + 8.0 / math.sqrt(df) + 0.5
    n_seg = 16
    total = 0.0
    step = hi / n_seg
    for i in range(n_seg):
        a = i * step
        b = a + step
        total += _integrate(integrand, a, b)
    return min(1.0, max(0.0, total))


def tukey_hsd(groups, labels=None):
    """Tukey HSD 多重比较（Tukey-Kramer，允许样本量不等）。
    groups: list[list[float]]，labels: 组名。
    返回 (anova_dict, comparisons)，comparisons 为
      list of dict(i, j, label_i, label_j, diff, q, p)。"""
    groups = [list(g) for g in groups]
    k = len(groups)
    if labels is None:
        labels = [str(i) for i in range(k)]
    av = oneway_anova(groups)
    msw = av['msw']
    df = av['df2']
    means = [mean(g) for g in groups]
    ns = [len(g) for g in groups]
    comps = []
    for i in range(k):
        for j in range(i + 1, k):
            diff = means[i] - means[j]
            se = math.sqrt(msw / 2.0 * (1.0 / ns[i] + 1.0 / ns[j])) if msw > 0 else 0.0
            q = abs(diff) / se if se > 0 else float('inf')
            p = 1.0 - ptukey(q, k, df) if se > 0 else 0.0
            p = min(1.0, max(0.0, p))
            comps.append({'i': i, 'j': j, 'label_i': labels[i], 'label_j': labels[j],
                          'diff': diff, 'q': q, 'p': p})
    return av, comps


def _p_max_abs_t(c, lams, df):
    """Dunnett 用：P(max_j |T_j| <= c)，T_j 为相关系数 lam_i*lam_j 的多元 t，自由度 df。
    先对公共正态分量 z 积分，再对 s ~ sqrt(chi2_df/df) 积分。"""
    if c <= 0:
        return 0.0
    rs = [math.sqrt(1.0 - l * l) for l in lams]

    def given_x(x):
        def fz(z):
            prod = _phi(z)
            for l, r in zip(lams, rs):
                prod *= _Phi((x - l * z) / r) - _Phi((-x - l * z) / r)
                if prod <= 0:
                    return 0.0
            return prod
        total = 0.0
        for i in range(4):
            a = -8.0 + 4.0 * i
            total += _integrate(fz, a, a + 4.0)
        return total

    if df > 2000:
        return min(1.0, max(0.0, given_x(c)))
    half_df = df / 2.0
    log_const = math.log(2.0) + half_df * math.log(half_df) - math.lgamma(half_df)

    def fs(s):
        if s <= 0:
            return 0.0
        return math.exp(log_const + (df - 1) * math.log(s) - half_df * s * s) * given_x(c * s)

    hi = 1.0 + 8.0 / math.sqrt(df) + 0.5
    n_seg = 8
    step = hi / n_seg
    total = sum(_integrate(fs, i * step, (i + 1) * step) for i in range(n_seg))
    return min(1.0, max(0.0, total))


def dunnett(control, others):
    """Dunnett 多重比较：各处理组与对照组比较，双侧，合并方差来自全部组的 ANOVA。
    control: list[float]，others: list[list[float]]。
    返回 (anova_dict, comparisons)，comparisons 为 list of dict(j, diff, t, p)，j 为 others 下标。"""
    av = oneway_anova([control] + list(others))
    msw, df = av['msw'], av['df2']
    n0 = len(control)
    m0 = mean(control)
    lams = [math.sqrt(len(g) / (len(g) + n0)) for g in others]
    comps = []
    for j, g in enumerate(others):
        diff = mean(g) - m0
        se = math.sqrt(msw * (1.0 / len(g) + 1.0 / n0)) if msw > 0 else 0.0
        if se > 0:
            t = diff / se
            p = 1.0 - _p_max_abs_t(abs(t), lams, df)
        else:
            t = 0.0 if diff == 0 else math.copysign(float('inf'), diff)
            p = 1.0 if diff == 0 else 0.0
        comps.append({'j': j, 'diff': diff, 't': t, 'p': min(1.0, max(0.0, p))})
    return av, comps


def holm_adjust(ps):
    """Holm 逐步校正。NaN 原样保留，只在有效 p 值之间校正，返回同序列表。"""
    valid = sorted((p, i) for i, p in enumerate(ps) if p == p)
    m = len(valid)
    out = list(ps)
    running = 0.0
    for rank, (p, i) in enumerate(valid):
        running = max(running, min(1.0, (m - rank) * p))
        out[i] = running
    return out


def compact_letter_display(labels, sig_pairs, order_by=None):
    """字母标记法（Piepho 2004 insert-and-absorb）。
    labels: 组名列表。
    sig_pairs: set of frozenset({i,j}) 表示 i,j 显著不同（不能共享字母）。
    order_by: 可选，按该数值降序给字母（通常按均值），使 'a' 给最大组。
    返回 dict{label: 字母串}。"""
    n = len(labels)
    # 每个 column 是一个集合，含本字母覆盖的组索引
    # 初始：一个 column 含所有组
    cols = [set(range(n))]

    def connected(group_set, idx_pairs):
        return True

    for pair in sig_pairs:
        i, j = tuple(pair)
        new_cols = []
        for col in cols:
            if i in col and j in col:
                # 违反：需要拆分
                c1 = set(col)
                c1.discard(j)
                c2 = set(col)
                c2.discard(i)
                new_cols.append(c1)
                new_cols.append(c2)
            else:
                new_cols.append(col)
        # 吸收：删去被其它 column 真包含的 column
        cols = _absorb(new_cols)

    # 去重并删空
    cols = [c for c in cols if c]
    # 字母排序：按 column 内最高 order 值排序，使重要组拿前面的字母
    if order_by is not None:
        cols.sort(key=lambda c: -max(order_by[i] for i in c))
    else:
        cols.sort(key=lambda c: min(c))

    letters = {lab: '' for lab in labels}
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    for ci, col in enumerate(cols):
        ch = alphabet[ci] if ci < 26 else alphabet[ci % 26] * (ci // 26 + 1)
        for idx in col:
            letters[labels[idx]] += ch
    # 保证每组至少一个字母
    for lab in labels:
        if letters[lab] == '':
            letters[lab] = alphabet[0]
    return letters


def _absorb(cols):
    """删去被其它集合真包含（子集）的集合。"""
    keep = []
    for i, c in enumerate(cols):
        if not c:
            continue
        subset = False
        for j, d in enumerate(cols):
            if i == j or not d:
                continue
            if c <= d and (c != d or i > j):
                subset = True
                break
        if not subset:
            keep.append(c)
    # 去重
    uniq = []
    for c in keep:
        if c not in uniq:
            uniq.append(c)
    return uniq


def p_to_stars(p):
    """p 值转显著性星号。"""
    if p != p:  # NaN
        return 'ns'
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return 'ns'


# =============================================================================
# 2. CFX reader（手写 xlsx 解析，兼容 CFX Maestro bug）
# =============================================================================

def _col_to_idx(col):
    """'A'->0, 'B'->1, 'AA'->26 ..."""
    idx = 0
    for ch in col:
        idx = idx * 26 + (ord(ch) - ord('A') + 1)
    return idx - 1


def _open_xlsx(path):
    """打开 xlsx，返回 (zipfile, namemap, sst)。
    兼容 CFX Maestro：ZIP 内部用反斜杠路径 + content_types/sharedStrings 大小写错误。"""
    z = zipfile.ZipFile(path)
    # 规范化条目名：反斜杠 -> 正斜杠，建 lower->真实名 映射
    namemap = {}
    for n in z.namelist():
        namemap[n.replace('\\', '/').lower()] = n
    sst = []
    key = next((k for k in namemap if k.endswith('sharedstrings.xml')), None)
    if key is not None:
        root = ET.fromstring(z.read(namemap[key]))
        for si in root:
            if _local(si.tag) != 'si':
                continue
            # 只取 <t>，跳过拼音注音 <rPh> 里的 <t>
            txt = ''.join(t.text or '' for t in si.iter()
                          if _local(t.tag) == 't' and not _in_rph(si, t))
            sst.append(txt)
    return z, namemap, sst


def _in_rph(si, target):
    for rph in si:
        if _local(rph.tag) == 'rPh' and any(x is target for x in rph.iter()):
            return True
    return False


def _parse_sheet_xml(data, sst):
    """解析一个 worksheet XML，返回 list[list]（按列顺序，缺列补 None）。"""
    root = ET.fromstring(data)
    rows = []
    for row in root.iter():
        if _local(row.tag) != 'row':
            continue
        cells = {}
        maxc = -1
        for c in row:
            if _local(c.tag) != 'c':
                continue
            ref = c.get('r', '')
            mcol = re.match(r'[A-Z]+', ref)
            if not mcol:
                continue
            ci = _col_to_idx(mcol.group())
            t = c.get('t')
            v = None
            for ch in c:
                lt = _local(ch.tag)
                if lt == 'v':
                    v = ch.text
                elif lt == 'is':
                    v = ''.join(x.text or '' for x in ch.iter() if _local(x.tag) == 't')
            if v is None:
                continue
            if t == 's':
                try:
                    val = sst[int(v)]
                except (ValueError, IndexError):
                    val = v
            elif t in ('str', 'inlineStr', 'e'):
                val = v
            elif t == 'b':
                val = (v == '1')
            else:
                try:
                    val = float(v)
                except (ValueError, TypeError):
                    val = v
            cells[ci] = val
            if ci > maxc:
                maxc = ci
        rows.append([cells.get(i) for i in range(maxc + 1)])
    return rows


def read_xlsx_rows(path):
    """读取 xlsx 第一个工作表，返回 list[list]（按列顺序，缺列补 None）。"""
    z, namemap, sst = _open_xlsx(path)
    try:
        # 找第一个 worksheet（按文件名里的数字排序）
        sheet_keys = [k for k in namemap
                      if '/worksheets/' in k and k.endswith('.xml') and 'rels' not in k]
        sheet_keys.sort(key=lambda k: _natural_key(k))
        if not sheet_keys:
            raise ValueError('xlsx 中找不到工作表')
        return _parse_sheet_xml(z.read(namemap[sheet_keys[0]]), sst)
    finally:
        z.close()


def read_xlsx_sheets(path):
    """按工作表名读取全部 sheet，返回 dict{sheet名: rows}。
    通过 workbook.xml 与其 rels 定位，兼容 Excel / WPS 另存后的文件。"""
    z, namemap, sst = _open_xlsx(path)
    try:
        wb_key = next((k for k in namemap if k.endswith('xl/workbook.xml')), None)
        rels_key = next((k for k in namemap if k.endswith('xl/_rels/workbook.xml.rels')), None)
        if wb_key is None or rels_key is None:
            raise ValueError('xlsx 结构不完整（缺 workbook.xml）')
        rid_target = {}
        for rel in ET.fromstring(z.read(namemap[rels_key])):
            rid_target[rel.get('Id')] = rel.get('Target', '')
        out = {}
        for el in ET.fromstring(z.read(namemap[wb_key])).iter():
            if _local(el.tag) != 'sheet':
                continue
            rid = next((v for k, v in el.attrib.items() if _local(k) == 'id'), None)
            target = rid_target.get(rid, '').replace('\\', '/').lower().lstrip('/')
            if not target.startswith('xl/'):
                target = 'xl/' + target
            real = namemap.get(target)
            if real is None:
                continue
            out[el.get('name')] = _parse_sheet_xml(z.read(real), sst)
        return out
    finally:
        z.close()


def read_csv_rows(path):
    """读取 CFX 导出的 CSV（UTF-8 / GBK 自动识别），数值列转 float。"""
    import csv
    with open(path, 'rb') as fh:
        raw = fh.read()
    for enc in ('utf-8-sig', 'gbk', 'latin-1'):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    rows = []
    for r in csv.reader(io.StringIO(text)):
        out = []
        for x in r:
            x = x.strip()
            if x == '':
                out.append(None)
                continue
            try:
                out.append(float(x))
            except ValueError:
                out.append(x)
        rows.append(out)
    return rows


def _local(tag):
    return tag.split('}')[-1]


def _natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', s)]


# CFX Maestro 的 .pcrd 是固定密码加密的 ZIP，内层为单个 XML
PCRD_PASSWORD = b'SecureCompressDecompressKeyiQ5V4Files!!##$$'

_PCRD_SAMPLE_TYPE = {'wcSample': 'Unkn', 'wcStandard': 'Std', 'wcNTC': 'NTC',
                     'wcPosCtrl': 'PosCtrl', 'wcNegCtrl': 'NegCtrl'}


def _read_pcrd_xml(path):
    """读取 pcrd 内层 XML。新版为加密 ZIP 包单个 XML，旧版可能直接是 XML 文本。"""
    with open(path, 'rb') as fh:
        raw = fh.read()
    if raw[:2] != b'PK':
        try:
            return raw.decode('utf-8-sig')
        except UnicodeDecodeError:
            raise ValueError('无法识别的 pcrd 文件（既不是加密 ZIP 也不是 XML）')
    z = zipfile.ZipFile(io.BytesIO(raw))
    try:
        names = [n for n in z.namelist() if not n.endswith('/')]
        if not names:
            raise ValueError('pcrd 压缩包为空')
        data = None
        for pwd in (None, PCRD_PASSWORD):
            try:
                data = z.read(names[0]) if pwd is None else z.read(names[0], pwd=pwd)
                break
            except RuntimeError:
                continue
        if data is None:
            raise ValueError('pcrd 解密失败（文件损坏或密码已变更）')
        return data.decode('utf-8-sig', 'replace')
    finally:
        z.close()


def parse_pcrd(path):
    """解析 CFX Maestro 原始 .pcrd 数据文件。
    孔位索引为行优先（0=A01, 1=A02, …）；Cq 取软件当前选中分析
    （dataAnalysisParameters 记录的孔组与步骤）下的 thresholdCycle。
    返回与 parse_cfx 相同的 list[dict(well, fluor, target, content, sample, cq)]。"""
    try:
        root = ET.fromstring(_read_pcrd_xml(path))
    except ET.ParseError as e:
        raise ValueError(f'pcrd 内容不是有效 XML：{e}')
    layers = list(root.iter('dyeLayer'))
    if not layers:
        raise ValueError('pcrd 中找不到板信息（dyeLayer）')
    fluors = [l.get('plateName') or '' for l in layers]
    if len([f for f in fluors if f]) > 1:
        raise ValueError(f'该 pcrd 含多个荧光通道（{"、".join(fluors)}），暂不支持直接导入；'
                         '请先在 CFX Maestro 中导出 Quantification Summary（xlsx/csv）再导入')
    rows, cols = int(layers[0].get('RowsCount') or 8), int(layers[0].get('ColumnsCount') or 12)
    # 孔注释：plateIndex -> 基因 / 样本 / 孔类型
    info = {}
    for wsel in root.iter('wellSample'):
        try:
            info[int(wsel.get('plateIndex'))] = wsel.attrib
        except (TypeError, ValueError):
            continue
    if not info:
        raise ValueError('pcrd 中找不到孔注释（wellSample）')
    # Cq：优先选「当前选中的孔组 + 选中的步骤」，退化为 All Wells 组或第一个有数据的分析
    dap = next(root.iter('dataAnalysisParameters'), None)
    sel_step = dap.get('selectedStepNumber') if dap is not None else None
    sel_group = dap.get('selectedWellGroupName') if dap is not None else None
    best = None
    for p in root.iter('dataAnalysisParam'):
        if next(p.iter('computedWellData'), None) is None:
            continue
        score = (p.get('WellGroupGUID') == sel_group, p.get('StepNumber') == sel_step,
                 p.get('wellGroupName') == 'All Wells')
        if best is None or score > best[0]:
            best = (score, p)
    cq_map = {}
    if best is not None:
        for cwd in best[1].iter('computedWellData'):
            if cwd.get('isComputed') != 'True':
                continue
            try:
                v = float(cwd.get('thresholdCycle'))
                idx = int(cwd.get('pIndex'))
            except (TypeError, ValueError):
                continue
            if math.isfinite(v):
                cq_map[idx] = v
    wells = []
    for idx in sorted(info):
        a = info[idx]
        stype = a.get('wellSampleType', '')
        if stype == 'wcEmpty':
            continue
        target = (a.get('geneName') or '').strip()
        sample = (a.get('conditionName') or a.get('sampleId') or '').strip()
        cq = cq_map.get(idx)
        if not target and not sample and cq is None:
            continue
        r, c = divmod(idx, cols)
        wells.append({'well': f'{chr(65 + r)}{c + 1:02d}',
                      'fluor': '/'.join(f for f in fluors if f),
                      'target': target,
                      'content': _PCRD_SAMPLE_TYPE.get(stype, stype),
                      'sample': sample,
                      'cq': cq})
    if not wells:
        raise ValueError('pcrd 中没有有效孔')
    return wells


def parse_cfx(path):
    """解析 CFX 文件（Quantification Summary / Cq Results 导出，或原始 .pcrd）。
    支持 .xlsx、.csv 与 .pcrd 三种格式。
    返回 list[dict(well, fluor, target, content, sample, cq)]，仅保留有效 Cq 的孔。"""
    low = str(path).lower()
    if low.endswith('.pcrd'):
        return parse_pcrd(path)
    rows = read_csv_rows(path) if low.endswith('.csv') else read_xlsx_rows(path)
    # 找表头行
    header_idx = None
    colmap = {}
    for i, r in enumerate(rows[:30]):
        vals = [str(x).strip() if x is not None else '' for x in r]
        if 'Cq' in vals and 'Sample' in vals and 'Target' in vals:
            header_idx = i
            for ci, name in enumerate(vals):
                if name:
                    colmap[name] = ci
            break
    if header_idx is None:
        raise ValueError('未找到 CFX 表头（需含 Target / Sample / Cq 列）')

    def cell(r, name):
        ci = colmap.get(name)
        if ci is None or ci >= len(r):
            return None
        return r[ci]

    wells = []
    for r in rows[header_idx + 1:]:
        if not r or all(x is None for x in r):
            continue
        cq = cell(r, 'Cq')
        target = cell(r, 'Target')
        sample = cell(r, 'Sample')
        # Cq 可能为空（无扩增 / NaN）
        cq_val = None
        if isinstance(cq, (int, float)):
            cq_val = float(cq)
        elif isinstance(cq, str):
            try:
                cq_val = float(cq)
            except ValueError:
                cq_val = None
        # 'NaN'/'inf' 文本或非有限值不是有效 Cq，按缺失处理，避免污染统计
        if cq_val is not None and not math.isfinite(cq_val):
            cq_val = None
        if target is None and sample is None and cq_val is None:
            continue
        wells.append({
            'well': cell(r, 'Well'),
            'fluor': cell(r, 'Fluor'),
            'target': str(target).strip() if target is not None else '',
            'content': cell(r, 'Content'),
            'sample': str(sample).strip() if sample is not None else '',
            'cq': cq_val,
        })
    return wells


# =============================================================================
# 3. 数据处理：技术重复分组 / QC / ΔΔCt
# =============================================================================

class ReplicateGroup:
    """一个样本×基因的技术重复组。"""

    def __init__(self, sample, target, cqs, wells=None):
        self.sample = sample
        self.target = target
        self.cqs_raw = [c for c in cqs if c is not None]
        self.wells = wells or []
        self.cqs_used = list(self.cqs_raw)
        self.outliers = []          # 被剔除的 Cq 值
        self.qc_pass = None
        self.qc_reason = ''

    @property
    def n_raw(self):
        return len(self.cqs_raw)

    @property
    def n_used(self):
        return len(self.cqs_used)

    @property
    def mean(self):
        return mean(self.cqs_used) if self.cqs_used else float('nan')

    @property
    def sd(self):
        return stdev(self.cqs_used) if self.n_used >= 2 else 0.0

    @property
    def cq_range(self):
        return (max(self.cqs_used) - min(self.cqs_used)) if self.cqs_used else 0.0


def group_replicates(wells):
    """把孔按 (sample, target) 分组成技术重复。返回 list[ReplicateGroup]。
    保持首次出现顺序。"""
    order = []
    buckets = {}
    for w in wells:
        if not w['target'] and not w['sample']:
            continue
        if w['cq'] is None:
            # 无 Cq 的孔（NTC 阴性 / 无扩增）仍记录但不参与
            pass
        key = (w['sample'], w['target'])
        if key not in buckets:
            buckets[key] = {'cqs': [], 'wells': []}
            order.append(key)
        if w['cq'] is not None:
            buckets[key]['cqs'].append(w['cq'])
        buckets[key]['wells'].append(w.get('well'))
    groups = []
    for key in order:
        sample, target = key
        b = buckets[key]
        groups.append(ReplicateGroup(sample, target, b['cqs'], b['wells']))
    return groups


def run_qc(groups, sd_threshold=0.5, metric='sd', remove_outlier=True):
    """技术重复合理性 QC。
    metric: 'sd'（标准差）或 'range'（极差 max-min）。
    remove_outlier=True 时，对 n>=3 且超阈的组尝试三选二剔除离群孔。
    原地修改 groups 的 qc 字段。返回统计 dict。"""
    n_pass = n_fail = n_fixed = 0
    for g in groups:
        if g.n_raw < 2:
            g.qc_pass = True
            g.qc_reason = '单孔，无法评估重复'
            n_pass += 1
            continue
        spread = _spread(g.cqs_used, metric)
        if spread <= sd_threshold:
            g.qc_pass = True
            g.qc_reason = 'OK'
            n_pass += 1
            continue
        # 超阈，尝试剔除离群
        if remove_outlier and g.n_raw >= 3:
            best = _best_drop_one(g.cqs_raw, metric)
            if best is not None and best['spread'] <= sd_threshold:
                g.cqs_used = best['kept']
                g.outliers = best['dropped']
                g.qc_pass = True
                g.qc_reason = '剔除离群孔后达标'
                n_fixed += 1
                n_pass += 1
                continue
        g.qc_pass = False
        g.qc_reason = '重复离散超阈'
        n_fail += 1
    return {'pass': n_pass, 'fail': n_fail, 'fixed': n_fixed, 'total': len(groups)}


def _spread(xs, metric):
    if len(xs) < 2:
        return 0.0
    if metric == 'range':
        return max(xs) - min(xs)
    return stdev(xs)


def _best_drop_one(cqs, metric):
    """从 cqs 中剔除 1 个孔，使剩余离散最小。返回最优方案。"""
    if len(cqs) < 3:
        return None
    best = None
    for i in range(len(cqs)):
        kept = cqs[:i] + cqs[i + 1:]
        sp = _spread(kept, metric)
        if best is None or sp < best['spread']:
            best = {'kept': kept, 'dropped': [cqs[i]], 'spread': sp}
    return best


def diagnose_dataset(wells):
    """加载数据后的概览体检，供界面提示用。返回各项计数与缺信息线索。"""
    targets = sorted({w['target'] for w in wells if w['target']})
    samples = sorted({w['sample'] for w in wells if w['sample']})
    grp_samples = {}
    for w in wells:
        if w['sample']:
            grp_samples.setdefault(auto_group_name(w['sample']), set()).add(w['sample'])
    groups = sorted(grp_samples)
    empty_cq = sum(1 for w in wells if w['cq'] is None)
    # 每个 (sample, target) 的孔数，找单孔（无技术重复）
    rep_count = {}
    for w in wells:
        if w['cq'] is not None and w['sample'] and w['target']:
            rep_count[(w['sample'], w['target'])] = rep_count.get((w['sample'], w['target']), 0) + 1
    single_well = sum(1 for n in rep_count.values() if n == 1)
    return {
        'n_wells': len(wells),
        'n_valid_cq': sum(1 for w in wells if w['cq'] is not None),
        'empty_cq': empty_cq,
        'targets': targets,
        'samples': samples,
        'groups': groups,
        'group_sizes': {g: len(s) for g, s in grp_samples.items()},
        'single_well_pairs': single_well,
    }


def auto_group_name(sample):
    """从样本名推断生物学分组名：仅当末尾数字前有明确分隔符(空格/-/_)时，
    才把数字当作重复编号剥离。
    'siNC-1'->'siNC'，'CK 2'->'CK'，'Liver-10'->'Liver'；
    不剥小数点(4.30->4.30)、不剥无分隔尾数(siRNA2->siRNA2)，避免误并不同条件。
    命名特殊时可在界面「自定义样本分组」手动指定。"""
    s = str(sample).strip()
    m = re.match(r'^(.*[^\s\-_.])[\s\-_]+\d+$', s)
    if m and m.group(1):
        return m.group(1).rstrip(' -_')
    return s


def compute_ddct(groups, ref_target, control_group, group_of=None,
                 sd_threshold=0.5, only_qc_pass=False):
    """ΔΔCt 相对定量。
    groups: QC 后的 ReplicateGroup 列表。
    ref_target: 内参基因名。
    control_group: 对照组名（calibrator）。
    group_of: dict{sample: group_name}，默认用 auto_group_name。
    返回 dict:
      targets        : 目的基因列表（不含内参）
      per_sample     : dict{(sample,target): {dct, ddct, rq, ref_cq, tgt_cq, group}}
      per_group      : dict{target: {group: {rqs:[...], mean, sem, n, dcts:[...]}}}
    """
    if group_of is None:
        group_of = {}
    # 建查表：sample -> {target: mean_cq}
    sample_target_cq = {}
    samples_order = []
    for g in groups:
        if only_qc_pass and not g.qc_pass:
            continue
        if g.n_used == 0:
            continue
        sample_target_cq.setdefault(g.sample, {})[g.target] = g.mean
        if g.sample not in samples_order:
            samples_order.append(g.sample)

    # 目的基因 = 除内参外所有 target
    all_targets = []
    for st in sample_target_cq.values():
        for t in st:
            if t != ref_target and t not in all_targets:
                all_targets.append(t)

    def grp(sample):
        return group_of.get(sample) or auto_group_name(sample)

    per_sample = {}
    skipped_no_ref = []        # 缺内参 Cq 而被跳过的样本
    skipped_no_target = {}     # target -> 缺该目的基因 Cq 的样本列表
    # 1) 每样本 ΔCt
    for sample in samples_order:
        st = sample_target_cq[sample]
        ref_cq = st.get(ref_target)
        if ref_cq is None:
            skipped_no_ref.append(sample)
            continue
        for target in all_targets:
            tgt_cq = st.get(target)
            if tgt_cq is None:
                skipped_no_target.setdefault(target, []).append(sample)
                continue
            dct = tgt_cq - ref_cq
            per_sample[(sample, target)] = {
                'sample': sample, 'target': target, 'group': grp(sample),
                'ref_cq': ref_cq, 'tgt_cq': tgt_cq, 'dct': dct,
                'ddct': None, 'rq': None,
            }

    # 2) 每基因对照组 ΔCt 均值作 calibrator
    calibrator = {}
    for target in all_targets:
        ctrl_dcts = [v['dct'] for v in per_sample.values()
                     if v['target'] == target and v['group'] == control_group]
        calibrator[target] = mean(ctrl_dcts) if ctrl_dcts else None

    # 3) ΔΔCt 与 RQ
    for key, v in per_sample.items():
        cal = calibrator[v['target']]
        if cal is None:
            continue
        v['ddct'] = v['dct'] - cal
        v['rq'] = 2.0 ** (-v['ddct'])

    # 4) 按基因×组汇总
    per_group = {}
    for target in all_targets:
        per_group[target] = {}
        groups_seen = []
        for v in per_sample.values():
            if v['target'] != target or v['rq'] is None:
                continue
            gname = v['group']
            if gname not in per_group[target]:
                per_group[target][gname] = {'rqs': [], 'dcts': [], 'samples': []}
                groups_seen.append(gname)
            per_group[target][gname]['rqs'].append(v['rq'])
            per_group[target][gname]['dcts'].append(v['dct'])
            per_group[target][gname]['samples'].append(v['sample'])
        for gname, d in per_group[target].items():
            d['mean'] = mean(d['rqs']) if d['rqs'] else float('nan')
            d['sem'] = sem(d['rqs']) if len(d['rqs']) >= 2 else 0.0
            d['sd'] = stdev(d['rqs']) if len(d['rqs']) >= 2 else 0.0
            d['n'] = len(d['rqs'])

    # 对照组在各目的基因上是否有数据（calibrator 为 None 即无效）
    control_missing = [t for t in all_targets if calibrator.get(t) is None]

    return {
        'targets': all_targets,
        'ref_target': ref_target,
        'control_group': control_group,
        'per_sample': per_sample,
        'per_group': per_group,
        'calibrator': calibrator,
        'skipped_no_ref': skipped_no_ref,
        'skipped_no_target': skipped_no_target,
        'control_missing': control_missing,
    }


METHOD_NAMES = {
    't-test': 'Student t 检验',
    'anova': '单因素 ANOVA + Tukey HSD',
    'dunnett': '单因素 ANOVA + Dunnett（各组 vs 对照）',
}


def run_stats(per_group_target, value_key='dcts', control_group=None, alpha=0.05, comp='tukey'):
    """对一个基因的各组做统计检验。均为参数检验，假定各组方差相等。
    per_group_target: dict{group: {dcts:[...], rqs:[...], ...}}。
    value_key: 用于检验的值（默认 'dcts'，统计学上更规范；也可用 'rqs'）。
    comp: ≥3 组时的多重比较方式，'tukey' 全部两两比较，'dunnett' 只比各组与对照
          （对照组不足 2 个重复时退回 Tukey）。
    返回 dict:
      method  : 't-test' / 'anova' / 'dunnett'
      n_groups
      p       : 总体 p（t 检验的 p 或 ANOVA 的 p）
      pairwise: list of dict(group_i, group_j, p, stars)   （Tukey 为全部两两，Dunnett 只含 vs 对照）
      letters : dict{group: 字母}  （Dunnett 时为空，字母法只适用全部两两比较）
      stars_vs_control: dict{group: stars}  （有对照组时，各组 vs 对照）
    """
    gnames = list(per_group_target.keys())
    data = [per_group_target[g][value_key] for g in gnames]
    # 过滤掉样本量不足的组
    valid = [(g, d) for g, d in zip(gnames, data) if len(d) >= 2]
    if len(valid) < 2:
        return {'method': None, 'n_groups': len(valid), 'p': float('nan'),
                'pairwise': [], 'letters': {}, 'stars_vs_control': {}}
    gnames = [v[0] for v in valid]
    data = [v[1] for v in valid]
    # 字母方向始终按图上展示量 RQ 排序，使 'a' 落到最高柱（与检验所用 value_key 无关）
    means_for_order = {}
    for g in gnames:
        rqs = per_group_target[g].get('rqs')
        means_for_order[g] = mean(rqs) if rqs else mean(per_group_target[g][value_key])

    result = {'n_groups': len(gnames), 'pairwise': [], 'letters': {},
              'stars_vs_control': {}}

    if len(gnames) == 2:
        t, df, p = student_ttest(data[0], data[1])
        result['method'] = 't-test'
        result['p'] = p
        result['pairwise'] = [{'group_i': gnames[0], 'group_j': gnames[1],
                               'p': p, 'stars': p_to_stars(p)}]
        # 字母（2 组也可给）
        sig = set()
        if p < alpha:
            sig.add(frozenset({0, 1}))
        result['letters'] = compact_letter_display(gnames, sig, order_by=[means_for_order[g] for g in gnames])
    elif comp == 'dunnett' and control_group in gnames:
        ci = gnames.index(control_group)
        others = [g for g in gnames if g != control_group]
        av, comps = dunnett(data[ci], [data[gnames.index(g)] for g in others])
        result['method'] = 'dunnett'
        result['p'] = av['p']
        result['anova'] = av
        for c in comps:
            result['pairwise'].append({'group_i': others[c['j']], 'group_j': control_group,
                                       'p': c['p'], 'stars': p_to_stars(c['p'])})
    else:
        av, comps = tukey_hsd(data, labels=gnames)
        result['method'] = 'anova'
        result['p'] = av['p']
        result['anova'] = av
        sig = set()
        idx = {g: i for i, g in enumerate(gnames)}
        for c in comps:
            c['stars'] = p_to_stars(c['p'])
            result['pairwise'].append({'group_i': c['label_i'], 'group_j': c['label_j'],
                                       'p': c['p'], 'stars': c['stars']})
            if c['p'] < alpha:
                sig.add(frozenset({idx[c['label_i']], idx[c['label_j']]}))
        result['letters'] = compact_letter_display(
            gnames, sig, order_by=[means_for_order[g] for g in gnames])

    # vs 对照组的星号：直接取主检验的两两比较结果，与字母法同源、不再自相矛盾
    if control_group and control_group in gnames:
        for g in gnames:
            if g == control_group:
                result['stars_vs_control'][g] = ''
                continue
            pc = next((c for c in result['pairwise']
                       if {c['group_i'], c['group_j']} == {g, control_group}), None)
            result['stars_vs_control'][g] = pc['stars'] if pc else 'ns'
    return result


# =============================================================================
# 4. xlsx writer（手写，带颜色标记）
# =============================================================================

# 标准样式索引：见 _STYLES_XML
S_DEFAULT = 0
S_HEADER = 1     # 加粗灰底
S_PASS = 2       # 绿底
S_FAIL = 3       # 红底
S_WARN = 4       # 黄底
S_NUM2 = 5       # 两位小数
S_NUM4 = 6       # 四位小数


def _cell_ref(col, row):
    s = ''
    col += 1
    while col > 0:
        col, rem = divmod(col - 1, 26)
        s = chr(ord('A') + rem) + s
    return f'{s}{row + 1}'


def write_xlsx(path, sheets):
    """写 xlsx。sheets: list of (name, rows)。
    rows: list[list[cell]]，cell 可为：
        值（str/数字/None） 或 (值, style_index) 元组。"""
    # 收集 shared strings
    sst = []
    sst_idx = {}

    def sid(s):
        if s not in sst_idx:
            sst_idx[s] = len(sst)
            sst.append(s)
        return sst_idx[s]

    sheet_xmls = []
    for sheet in sheets:
        name, rows = sheet[0], sheet[1]
        opts = sheet[2] if len(sheet) > 2 else {}
        buf = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
               '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">']
        if opts.get('freeze_header'):
            buf.append('<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" '
                       'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>')
        if opts.get('col_widths'):
            buf.append('<cols>')
            for ci, w in enumerate(opts['col_widths']):
                buf.append(f'<col min="{ci + 1}" max="{ci + 1}" width="{w}" customWidth="1"/>')
            buf.append('</cols>')
        buf.append('<sheetData>')
        for ri, row in enumerate(rows):
            buf.append(f'<row r="{ri + 1}">')
            for ci, cell in enumerate(row):
                style = 0
                val = cell
                if isinstance(cell, tuple):
                    val, style = cell
                ref = _cell_ref(ci, ri)
                if val is None or val == '':
                    if style:
                        buf.append(f'<c r="{ref}" s="{style}"/>')
                    continue
                if isinstance(val, bool):
                    val = str(val)
                if isinstance(val, (int, float)) and not isinstance(val, bool):
                    if val != val or val in (float('inf'), float('-inf')):
                        buf.append(f'<c r="{ref}" s="{style}" t="inlineStr"><is><t>{_xml_text(val)}</t></is></c>')
                    else:
                        buf.append(f'<c r="{ref}" s="{style}"><v>{val}</v></c>')
                else:
                    i = sid(str(val))
                    buf.append(f'<c r="{ref}" s="{style}" t="s"><v>{i}</v></c>')
            buf.append('</row>')
        buf.append('</sheetData>')
        if opts.get('autofilter') and rows and rows[0]:
            buf.append(f'<autoFilter ref="A1:{_cell_ref(len(rows[0]) - 1, max(len(rows) - 1, 0))}"/>')
        buf.append('</worksheet>')
        sheet_xmls.append(''.join(buf))

    # sharedStrings.xml
    sst_xml = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
               f'<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" count="{len(sst)}" uniqueCount="{len(sst)}">']
    for s in sst:
        sst_xml.append(f'<si><t xml:space="preserve">{_xml_text(s)}</t></si>')
    sst_xml.append('</sst>')
    sst_xml = ''.join(sst_xml)

    # workbook.xml
    wb = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
          '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
          'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>']
    for i, sheet in enumerate(sheets):
        safe = _xml_escape(sheet[0][:31])
        wb.append(f'<sheet name="{safe}" sheetId="{i + 1}" r:id="rId{i + 1}"/>')
    wb.append('</sheets></workbook>')
    wb = ''.join(wb)

    # workbook rels
    rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
    for i in range(len(sheets)):
        rels.append(f'<Relationship Id="rId{i + 1}" '
                    f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
                    f'Target="worksheets/sheet{i + 1}.xml"/>')
    n = len(sheets)
    rels.append(f'<Relationship Id="rId{n + 1}" '
                f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/sharedStrings" '
                f'Target="sharedStrings.xml"/>')
    rels.append(f'<Relationship Id="rId{n + 2}" '
                f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
                f'Target="styles.xml"/>')
    rels.append('</Relationships>')
    rels = ''.join(rels)

    # content types
    ct = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
          '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
          '<Default Extension="xml" ContentType="application/xml"/>',
          '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
          '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
          '<Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/>']
    for i in range(len(sheets)):
        ct.append(f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" '
                  f'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
    ct.append('</Types>')
    ct = ''.join(ct)

    root_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                 '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                 '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
                 '</Relationships>')

    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml', ct)
        z.writestr('_rels/.rels', root_rels)
        z.writestr('xl/workbook.xml', wb)
        z.writestr('xl/_rels/workbook.xml.rels', rels)
        z.writestr('xl/sharedStrings.xml', sst_xml)
        z.writestr('xl/styles.xml', _STYLES_XML)
        for i, sx in enumerate(sheet_xmls):
            z.writestr(f'xl/worksheets/sheet{i + 1}.xml', sx)


# styles.xml：定义编号格式、字体、填充色、对应 cellXfs
_STYLES_XML = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<numFmts count="2">
<numFmt numFmtId="164" formatCode="0.00"/>
<numFmt numFmtId="165" formatCode="0.0000"/>
</numFmts>
<fonts count="2">
<font><sz val="11"/><name val="Calibri"/></font>
<font><b/><sz val="11"/><name val="Calibri"/></font>
</fonts>
<fills count="5">
<fill><patternFill patternType="none"/></fill>
<fill><patternFill patternType="gray125"/></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FFD9EAD3"/></patternFill></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FFF4CCCC"/></patternFill></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FFFFF2CC"/></patternFill></fill>
</fills>
<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="7">
<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
<xf numFmtId="0" fontId="1" fillId="1" borderId="0" xfId="0" applyFont="1" applyFill="1"/>
<xf numFmtId="0" fontId="0" fillId="2" borderId="0" xfId="0" applyFill="1"/>
<xf numFmtId="0" fontId="0" fillId="3" borderId="0" xfId="0" applyFill="1"/>
<xf numFmtId="0" fontId="0" fillId="4" borderId="0" xfId="0" applyFill="1"/>
<xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>
<xf numFmtId="165" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>
</cellXfs>
<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>'''


# =============================================================================
# 5. SVG 绘图（期刊级柱状图）
# =============================================================================

# 科研友好配色 —— Nature/Cell 极简风：低饱和、协调、克制
# 柔和谱系内的低饱和色，避免高饱和原色，打印与色弱友好
PALETTE = ['#6E8CAD', '#D49A66', '#7FB3A3', '#BC7C8C', '#988BB5',
           '#AE9277', '#CC9BBD', '#9DA2A8', '#CBB682', '#86B2C2']

# 绘图视觉常量（集中管理，便于统一风格）
_INK = '#1A1A1A'          # 主文字/标注（近黑而非纯黑）
_INK_SOFT = '#5A5F66'     # 次要文字（刻度数字）
_AXIS = '#33373B'         # 轴线（深灰）
_GRID = '#EEEFF1'         # 极淡水平参考线
_ERR = '#3A3F44'          # 误差棒
_DOT = '#33383D'          # 散点描边/填充


def _hex_mix(hex_color, other, t):
    """在 hex_color 与 other 之间按 t 线性插值（t=0 取 hex_color）。返回 #rrggbb。"""
    def _c(h):
        h = h.lstrip('#')
        return (int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16))
    r1, g1, b1 = _c(hex_color)
    r2, g2, b2 = _c(other)
    r = round(r1 + (r2 - r1) * t)
    g = round(g1 + (g2 - g1) * t)
    b = round(b1 + (b2 - b1) * t)
    return '#%02X%02X%02X' % (r, g, b)


def make_barplot_geometry(groups, means, errs, points=None, letters=None,
                          stars=None, control_group=None,
                          width=520, height=400, title='', ylabel='Relative expression (2^-ΔΔCt)',
                          xlabel=''):
    """计算柱状图几何（坐标），供 SVG / Canvas 共享渲染。
    groups: 组名列表; means/errs: 等长数值; points: list[list[float]] 各组散点。
    返回 dict 描述所有图元。"""
    n = len(groups)
    # 期刊版式：左留白容下旋转的 y 轴标题，右侧极窄；有标题时顶部稍宽
    margin = {'left': 70, 'right': 26, 'top': 40 if title else 26, 'bottom': 52}
    plot_w = width - margin['left'] - margin['right']
    plot_h = height - margin['top'] - margin['bottom']

    # y 轴范围：先按数据极值定刻度，再让上界落在“漂亮刻度”上，留出标注空间
    all_top = []
    for i in range(n):
        all_top.append(means[i] + (errs[i] if errs else 0))
        if points and i < len(points):
            all_top.extend(points[i])
    data_max = max(all_top) if all_top else 1.0
    if data_max <= 0:
        data_max = 1.0
    ymin = 0.0
    # 顶部留 ~14% 给显著性字母/星号，再吸附到整齐刻度
    ticks = _nice_ticks(ymin, data_max * 1.14, 5)
    ymax = max(ticks[-1], data_max * 1.06)

    def y2px(v):
        return margin['top'] + plot_h * (1 - (v - ymin) / (ymax - ymin))

    band = plot_w / n
    bar_w = min(band * 0.46, 56)   # 略瘦的柱，留白更讲究
    bars = []
    for i in range(n):
        cx = margin['left'] + band * (i + 0.5)
        x0 = cx - bar_w / 2
        top = y2px(means[i])
        base = y2px(0)
        col = PALETTE[i % len(PALETTE)]
        pts = []
        if points and i < len(points):
            pys = [y2px(v) for v in points[i]]
            xs = _jitter_offsets(pys, bar_w * 0.66)
            pts = [(cx + dx, py) for dx, py in zip(xs, pys)]
        err_top = y2px(means[i] + errs[i]) if errs else top
        err_bot = y2px(max(0, means[i] - errs[i])) if errs else base
        # 标注基准取“柱顶/误差棒顶/最高散点”三者最高处（像素最小 y）
        cand = [err_top, top]
        if pts:
            cand.append(min(p[1] for p in pts))
        bars.append({
            'group': groups[i], 'cx': cx, 'x0': x0, 'bar_w': bar_w,
            'top': top, 'base': base, 'color': col,
            'mean': means[i], 'err': errs[i] if errs else 0,
            'err_top': err_top, 'err_bot': err_bot,
            'points': pts,
            'letter': letters.get(groups[i]) if letters else None,
            'star': stars.get(groups[i]) if stars else None,
            'label_y': min(cand),
        })
    return {
        'width': width, 'height': height, 'margin': margin,
        'plot_w': plot_w, 'plot_h': plot_h,
        'ymin': ymin, 'ymax': ymax, 'ticks': ticks, 'y2px': y2px,
        'bars': bars, 'title': title, 'ylabel': ylabel, 'xlabel': xlabel,
    }


def _jitter_offsets(pys, span):
    """给定散点的纵坐标(像素)，返回每点的横向偏移，使纵向接近的点左右错开避免重叠。
    确定性、对称分布；不依赖随机数，渲染可复现。span 为允许的总横向宽度。"""
    n = len(pys)
    if n == 0:
        return []
    if n == 1:
        return [0.0]
    # 按 y 排序后分箱：纵向距离 < r_collide 视为同层，同层内左右对称展开
    r_collide = 5.4
    order = sorted(range(n), key=lambda k: pys[k])
    offs = [0.0] * n
    layer = []

    def flush(layer_idxs):
        m = len(layer_idxs)
        if m == 1:
            offs[layer_idxs[0]] = 0.0
            return
        # 对称：-(m-1)/2 .. (m-1)/2，步距受 span 约束
        step = min(span / max(m - 1, 1), 6.4)
        start = -step * (m - 1) / 2.0
        for k, gi in enumerate(layer_idxs):
            offs[gi] = start + step * k

    prev_y = None
    for gi in order:
        y = pys[gi]
        if prev_y is None or (y - prev_y) <= r_collide:
            layer.append(gi)
        else:
            flush(layer)
            layer = [gi]
        prev_y = y
    if layer:
        flush(layer)
    return offs


def _nice_ticks(lo, hi, n=5):
    """生成漂亮的刻度值。"""
    rng = hi - lo
    if rng <= 0:
        return [0, 1]
    raw = rng / n
    mag = 10 ** math.floor(math.log10(raw))
    norm = raw / mag
    if norm < 1.5:
        step = 1
    elif norm < 3:
        step = 2
    elif norm < 7:
        step = 5
    else:
        step = 10
    step *= mag
    ticks = []
    v = 0.0
    while v <= hi + step * 0.5:
        if v >= lo - step * 0.5:
            ticks.append(round(v, 10))
        v += step
    return ticks


def geometry_to_svg(geo, font='Arial, Helvetica, sans-serif'):
    """把几何渲染为 SVG 字符串（Nature/Cell 极简风：细轴、克制配色、精确标注）。"""
    W, H = geo['width'], geo['height']
    m = geo['margin']
    y2px = geo['y2px']
    base = y2px(0)
    s = []
    s.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
             f'viewBox="0 0 {W} {H}" font-family="{font}">')
    s.append(f'<rect width="{W}" height="{H}" fill="#FFFFFF"/>')

    # ---- 标题：轻量字重、居中、近黑，不喧宾夺主 ----
    if geo['title']:
        s.append(f'<text x="{W/2:.1f}" y="22" font-size="14" font-weight="500" '
                 f'letter-spacing="0.2" text-anchor="middle" fill="{_INK}">'
                 f'{_xml_escape(geo["title"])}</text>')

    # ---- 极淡水平参考线（仅在刻度处，>0；不画到 0 线以免与 x 轴重叠）----
    for t in geo['ticks']:
        if t <= geo['ymin']:
            continue
        yp = y2px(t)
        s.append(f'<line x1="{m["left"]:.1f}" y1="{yp:.1f}" x2="{W-m["right"]:.1f}" y2="{yp:.1f}" '
                 f'stroke="{_GRID}" stroke-width="1"/>')

    # ---- y 轴刻度文字 + 外向短刻度线 ----
    for t in geo['ticks']:
        yp = y2px(t)
        s.append(f'<text x="{m["left"]-9:.1f}" y="{yp+3.6:.1f}" font-size="11" '
                 f'text-anchor="end" fill="{_INK_SOFT}">{_fmt_tick(t)}</text>')
        s.append(f'<line x1="{m["left"]-4:.1f}" y1="{yp:.1f}" x2="{m["left"]:.1f}" y2="{yp:.1f}" '
                 f'stroke="{_AXIS}" stroke-width="1"/>')

    # ---- 坐标轴：细而清晰，仅左/下两条 ----
    s.append(f'<line x1="{m["left"]:.1f}" y1="{m["top"]:.1f}" x2="{m["left"]:.1f}" y2="{base:.1f}" '
             f'stroke="{_AXIS}" stroke-width="1"/>')
    s.append(f'<line x1="{m["left"]:.1f}" y1="{base:.1f}" x2="{W-m["right"]:.1f}" y2="{base:.1f}" '
             f'stroke="{_AXIS}" stroke-width="1"/>')

    # ---- y 轴标题：把 2^-ΔΔCt 渲染成真正的上标 ----
    # cairosvg 在 text-anchor=middle 下 tspan 上标会错位，故用 start 锚点手动居中：
    # 锚点放在左缘 x=yl_x、y=(cy+半宽)，绕锚点旋转 -90°，文本向上展开恰好居中于 cy
    cy = m['top'] + geo['plot_h'] / 2
    yl_fs = 12.5
    yl_x = 19
    yl_w = _approx_text_width(geo['ylabel'], yl_fs)
    yl_anchor_y = cy + yl_w / 2.0
    s.append(f'<text x="{yl_x}" y="{yl_anchor_y:.1f}" font-size="{yl_fs}" fill="{_INK}" '
             f'transform="rotate(-90 {yl_x} {yl_anchor_y:.1f})">'
             f'{_ylabel_tspans(geo["ylabel"])}</text>')

    # ---- 柱 + 误差棒 + 散点 + 标注 ----
    for b in geo['bars']:
        h = base - b['top']
        col = b['color']
        # 柔和实心柱 + 同色系更深的细描边，干净不刺眼
        edge = _hex_mix(col, _INK, 0.30)
        s.append(f'<rect x="{b["x0"]:.2f}" y="{b["top"]:.2f}" width="{b["bar_w"]:.2f}" '
                 f'height="{max(0,h):.2f}" fill="{col}" fill-opacity="0.92" '
                 f'stroke="{edge}" stroke-width="0.8" rx="1"/>')

        # 误差棒：细、上须为主（SEM 通常单向更显克制），深灰
        if b['err'] > 0:
            cx = b['cx']
            cap = min(b['bar_w'] * 0.20, 6.0)
            s.append(f'<line x1="{cx:.2f}" y1="{b["err_bot"]:.2f}" x2="{cx:.2f}" y2="{b["err_top"]:.2f}" '
                     f'stroke="{_ERR}" stroke-width="1"/>')
            s.append(f'<line x1="{cx-cap:.2f}" y1="{b["err_top"]:.2f}" x2="{cx+cap:.2f}" y2="{b["err_top"]:.2f}" '
                     f'stroke="{_ERR}" stroke-width="1"/>')
            s.append(f'<line x1="{cx-cap:.2f}" y1="{b["err_bot"]:.2f}" x2="{cx+cap:.2f}" y2="{b["err_bot"]:.2f}" '
                     f'stroke="{_ERR}" stroke-width="1"/>')

        # 散点：小而雅，白描边让重叠点可分辨
        for (px, py) in b['points']:
            s.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="2.1" '
                     f'fill="{_DOT}" fill-opacity="0.62" stroke="#FFFFFF" stroke-width="0.6"/>')

        # 显著性标注（字母优先，其次星号）：贴最高图元上方，位置精准
        lab = b['letter'] if b['letter'] else b['star']
        if lab:
            ly = b['label_y'] - 9
            is_letter = bool(b['letter']) and not b['star']
            fs = 12 if is_letter else 13
            fw = '600' if is_letter else '700'
            # 星号竖直居中略微下沉以贴近 baseline 视觉
            dy = 0 if is_letter else 1
            s.append(f'<text x="{b["cx"]:.2f}" y="{ly+dy:.2f}" font-size="{fs}" font-weight="{fw}" '
                     f'text-anchor="middle" fill="{_INK}">{_xml_escape(str(lab))}</text>')

        # x 轴组名
        s.append(f'<text x="{b["cx"]:.2f}" y="{base+16:.1f}" font-size="11.5" '
                 f'text-anchor="middle" fill="{_INK}">{_xml_escape(str(b["group"]))}</text>')

    if geo['xlabel']:
        s.append(f'<text x="{m["left"]+geo["plot_w"]/2:.1f}" y="{H-6}" font-size="12" '
                 f'text-anchor="middle" fill="{_INK}">{_xml_escape(geo["xlabel"])}</text>')
    # 角落署名（中文字体栈，Windows/Office 正常显示）
    s.append(f'<text x="{W-6}" y="{H-5}" font-size="8" text-anchor="end" '
             f'font-family="Microsoft YaHei, PingFang SC, SimHei, sans-serif" '
             f'fill="#BBBBBB">自动挡赛车手制作</text>')
    s.append('</svg>')
    return '\n'.join(s)


def _approx_text_width(label, font_size):
    """估算字符串在给定字号下的渲染宽度（用于 start 锚点手动居中）。
    Arial 经验平均字宽 ≈ 0.52em；上标 '2^...' 中的 ^ 不计宽，上标字小算 0.32em。"""
    w = 0.0
    sup = False
    i = 0
    while i < len(label):
        ch = label[i]
        if ch == '^':
            sup = True
            i += 1
            continue
        if ch in ')':
            sup = False
        em = 0.32 if sup else (0.30 if ch == ' ' else 0.52)
        # 宽字母略宽，窄字符略窄
        if ch in 'mwMW':
            em += 0.18
        elif ch in 'iljtIfr.,()':
            em -= 0.18
        w += em
        i += 1
    return w * font_size


def _ylabel_tspans(label):
    """把 y 轴标题里的 '2^-ΔΔCt' 渲染为真正上标；其余原样转义。
    识别形如 '...(2^-ΔΔCt)...' 段，无则整体原样输出。
    用 dy 偏移实现上标（cairosvg 对 baseline-shift 支持不稳）并在 start 锚点下显式复位。"""
    m = re.search(r'2\^(-?[^\s)]+)', label)
    if not m:
        return _xml_escape(label)
    pre = label[:m.start()]
    expo = m.group(1)
    tail = label[m.end():]
    return (f'{_xml_escape(pre)}2'
            f'<tspan dy="-4.5" font-size="8.5">{_xml_escape(expo)}</tspan>'
            f'<tspan dy="4.5">{_xml_escape(tail)}</tspan>')


def _fmt_tick(v):
    if v == int(v):
        return str(int(v))
    return f'{v:g}'
