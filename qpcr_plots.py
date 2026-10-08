# -*- coding: utf-8 -*-
"""
qPCR 数据中心 · 出图模板（纯标准库生成 SVG）

坐标单位为 pt（1 mm = 72/25.4 pt），SVG 宽高直接写成 mm，插入 Word / AI / PPT 即为实际尺寸。
模板：
  bar      柱状 + 个体点（均值 ± SEM/SD）
  box      箱线 + 个体点
  dot      均值点图（均值横线 + 误差须 + 个体点）
  grouped  多基因分组柱状图
  heatmap  log2FC 热图
  line     时间序列折线图
"""

import math
from xml.sax.saxutils import escape as _esc

import qpcr_core as q

PT_PER_MM = 72 / 25.4

# ---------------------------------------------------------------- 配色
PALETTES = {
    'nature':  ('Nature 极简', ['#6E8CAD', '#D49A66', '#7FB3A3', '#BC7C8C', '#988BB5', '#AE9277', '#CC9BBD', '#9DA2A8', '#CBB682', '#86B2C2']),
    'npg':     ('NPG 经典', ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85']),
    'aaas':    ('Science', ['#3B4992', '#EE0000', '#008B45', '#631879', '#008280', '#BB0021', '#5F559B', '#A20056', '#808180', '#1B1919']),
    'lancet':  ('Lancet', ['#00468B', '#ED0000', '#42B540', '#0099B4', '#925E9F', '#FDAF91', '#AD002A', '#ADB6B6', '#1B1919']),
    'nejm':    ('NEJM', ['#BC3C29', '#0072B5', '#E18727', '#20854E', '#7876B1', '#6F99AD', '#FFDC91', '#EE4C97']),
    'jco':     ('JCO', ['#0073C2', '#EFC000', '#868686', '#CD534C', '#7AA6DC', '#003C67', '#8F7700', '#3B3B3B']),
    'okabe':   ('色盲友好', ['#0072B2', '#E69F00', '#009E73', '#D55E00', '#56B4E9', '#CC79A7', '#F0E442', '#999999']),
    'morandi': ('莫兰迪', ['#8E9AAF', '#C9ADA7', '#9DB0A3', '#D8C3A5', '#9A8C98', '#C1A192', '#A3B1C2', '#BFB5AF']),
    'pastel':  ('日系淡彩', ['#A7C7E7', '#F4B6C2', '#B5DDC4', '#F8D49B', '#C9B6E4', '#F2C6A0', '#A8D8D8', '#D9D3C7']),
    'gray':    ('黑白灰', ['#3A3A3A', '#7A7A7A', '#AFAFAF', '#D6D6D6', '#545454', '#959595']),
}
HEATMAPS = {
    'rdbu':     ('蓝白红', ['#2166AC', '#F7F7F7', '#B2182B']),
    'coolwarm': ('冷暖', ['#3B4CC0', '#EDEDED', '#B40426']),
    'puor':     ('紫白橙', ['#542788', '#F7F7F7', '#B35806']),
    'brbg':     ('棕白绿', ['#8C510A', '#F5F5F5', '#01665E']),
    'pigr':     ('粉白绿', ['#C51B7D', '#F7F7F7', '#4D9221']),
}
PRESETS = [
    {'key': 'nature',  'label': 'Nature 极简', 'palette': 'nature',  'fill': 'solid',   'point': 'dark',  'err': 'sem', 'ctrl_gray': False, 'heat': 'rdbu'},
    {'key': 'cell',    'label': 'Cell 灰对照', 'palette': 'npg',     'fill': 'solid',   'point': 'dark',  'err': 'sem', 'ctrl_gray': True,  'heat': 'coolwarm'},
    {'key': 'science', 'label': 'Science 经典', 'palette': 'aaas',   'fill': 'light',   'point': 'color', 'err': 'sd',  'ctrl_gray': False, 'heat': 'rdbu'},
    {'key': 'lancet',  'label': 'Lancet 描边', 'palette': 'lancet',  'fill': 'outline', 'point': 'color', 'err': 'sd',  'ctrl_gray': False, 'heat': 'coolwarm'},
    {'key': 'okabe',   'label': '色盲友好',     'palette': 'okabe',   'fill': 'solid',   'point': 'dark',  'err': 'sem', 'ctrl_gray': False, 'heat': 'puor'},
    {'key': 'morandi', 'label': '莫兰迪',       'palette': 'morandi', 'fill': 'solid',   'point': 'open',  'err': 'sem', 'ctrl_gray': False, 'heat': 'brbg'},
    {'key': 'pastel',  'label': '日系淡彩',     'palette': 'pastel',  'fill': 'light',   'point': 'color', 'err': 'sem', 'ctrl_gray': True,  'heat': 'pigr'},
    {'key': 'mono',    'label': '黑白印刷',     'palette': 'gray',    'fill': 'outline', 'point': 'dark',  'err': 'sd',  'ctrl_gray': False, 'heat': 'rdbu'},
]
TEMPLATES = [
    ('bar', '柱状 + 散点'), ('box', '箱线 + 散点'), ('dot', '均值点图'),
    ('grouped', '多基因分组柱'), ('heatmap', 'log2FC 热图'), ('line', '时间序列折线'),
]
FONTS = {'arial': 'Arial, Helvetica, sans-serif', 'helvetica': 'Helvetica, Arial, sans-serif',
         'times': "'Times New Roman', Times, serif"}

DEFAULTS = {
    'template': 'bar', 'gene': '', 'genes': [], 'groups': [],
    'preset': 'nature', 'palette': 'nature', 'heat': 'rdbu', 'colors': {},
    'fill': 'solid', 'point': 'dark', 'err': 'sem', 'ctrl_gray': False,
    'y': 'rq', 'annot': 'auto', 'width_mm': 89, 'height_mm': 70, 'font_pt': 7,
    'font': 'arial', 'title': '', 'ylabel': '', 'xlabel': '', 'italic': True,
    'heat_values': True, 'signature': False, 'bar_width': 0.62,
}

INK = '#1A1A1A'
INK_SOFT = '#4A4A4A'
CTRL_GRAY = '#A6A6A6'
LW = 0.5           # 轴线 / 误差棒线宽 pt


# ---------------------------------------------------------------- 工具
def _hex(c):
    c = c.lstrip('#')
    return tuple(int(c[i:i + 2], 16) for i in (0, 2, 4))


def mix(a, b, t):
    ra, ga, ba = _hex(a)
    rb, gb, bb = _hex(b)
    return '#%02X%02X%02X' % (round(ra + (rb - ra) * t), round(ga + (gb - ga) * t), round(ba + (bb - ba) * t))


def ramp(stops, t):
    """三色发散色带取色，t ∈ [-1, 1]。"""
    t = max(-1.0, min(1.0, t))
    return mix(stops[1], stops[0], -t) if t < 0 else mix(stops[1], stops[2], t)


def text_width(s, fs):
    w = 0.0
    for ch in str(s):
        o = ord(ch)
        if o > 0x2E80:
            w += 1.0
        elif ch in 'ijlI.,:;|!\'':
            w += 0.26
        elif ch in 'ft()[]{} -/':
            w += 0.33
        elif ch in 'mwMW':
            w += 0.83
        elif ch.isupper():
            w += 0.67
        else:
            w += 0.55
    return w * fs


def _rich_parts(s):
    """解析 ^{上标} 与 _{下标} 标记，返回 [(文本, 'n'|'sup'|'sub')]。"""
    parts, i, buf = [], 0, ''
    s = str(s)
    while i < len(s):
        if s[i] in '^_' and i + 1 < len(s) and s[i + 1] == '{':
            j = s.find('}', i + 2)
            if j > 0:
                if buf:
                    parts.append((buf, 'n'))
                    buf = ''
                parts.append((s[i + 2:j], 'sup' if s[i] == '^' else 'sub'))
                i = j + 1
                continue
        buf += s[i]
        i += 1
    if buf:
        parts.append((buf, 'n'))
    return parts


def rich_width(s, fs):
    return sum(text_width(t, fs * (0.7 if k != 'n' else 1)) for t, k in _rich_parts(s))


def nice_step(raw):
    if raw <= 0:
        return 1.0
    mag = 10 ** math.floor(math.log10(raw))
    norm = raw / mag
    for m in (1, 2, 2.5, 5, 10):
        if norm <= m:
            return m * mag
    return 10 * mag


def nice_ticks(lo, hi, n=5):
    if hi <= lo:
        hi = lo + 1
    step = nice_step((hi - lo) / n)
    start = math.floor(lo / step + 1e-9) * step
    end = math.ceil(hi / step - 1e-9) * step
    ticks, v = [], start
    while v <= end + step * 0.01:
        ticks.append(round(v, 10))
        v += step
    return ticks


def fmt_tick(v):
    if abs(v - round(v)) < 1e-9:
        return str(int(round(v))).replace('-', '−')
    return f'{v:g}'.replace('-', '−')


def quantile(xs, p):
    xs = sorted(xs)
    if not xs:
        return float('nan')
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def jitter(pys, span, r):
    """纵向相近的点左右对称错开，确定性、可复现。"""
    n = len(pys)
    if n <= 1:
        return [0.0] * n
    order = sorted(range(n), key=lambda k: pys[k])
    offs = [0.0] * n
    layer, prev = [], None

    def flush(ids):
        m = len(ids)
        if m == 1:
            offs[ids[0]] = 0.0
            return
        step = min(span / (m - 1), r * 2.3)
        for k, gi in enumerate(ids):
            offs[gi] = -step * (m - 1) / 2 + step * k

    for gi in order:
        if prev is None or pys[gi] - prev <= r * 2.1:
            layer.append(gi)
        else:
            flush(layer)
            layer = [gi]
        prev = pys[gi]
    flush(layer)
    return offs


class Svg:
    def __init__(self, w, h, font, fs):
        self.w, self.h, self.font, self.fs = w, h, font, fs
        self.el = []

    def add(self, s):
        self.el.append(s)

    def line(self, x1, y1, x2, y2, stroke=INK, sw=LW, dash=None, cap='butt'):
        d = f' stroke-dasharray="{dash}"' if dash else ''
        self.add(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{stroke}" '
                 f'stroke-width="{sw}" stroke-linecap="{cap}"{d}/>')

    def rect(self, x, y, w, h, fill, stroke='none', sw=LW):
        if h < 0:
            y, h = y + h, -h
        self.add(f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(w, 0):.2f}" height="{max(h, 0):.2f}" '
                 f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')

    def circle(self, x, y, r, fill, stroke='none', sw=0.3, op=1.0):
        o = f' fill-opacity="{op}"' if op < 1 else ''
        self.add(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{fill}"{o} stroke="{stroke}" stroke-width="{sw}"/>')

    def path(self, d, stroke, sw=0.8, fill='none'):
        self.add(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" '
                 f'stroke-linejoin="round" stroke-linecap="round"/>')

    def text(self, x, y, s, size=None, anchor='middle', weight=None, italic=False,
             fill=INK, rotate=None):
        size = size or self.fs
        parts = _rich_parts(s)
        attrs = f' font-size="{size:.2f}" fill="{fill}"'
        if weight:
            attrs += f' font-weight="{weight}"'
        if italic:
            attrs += ' font-style="italic"'
        rich = any(k != 'n' for _, k in parts)
        if rich and anchor != 'start':
            # 富文本统一用 start 锚点手动居中，避免各渲染器对 middle + tspan 处理不一
            w = rich_width(s, size)
            if rotate is None:
                x = x - (w / 2 if anchor == 'middle' else w)
            else:
                y = y + (w / 2 if anchor == 'middle' else w)
            anchor = 'start'
        tr = f' transform="rotate({rotate} {x:.2f} {y:.2f})"' if rotate is not None else ''
        if not rich:
            body = _esc(str(s))
        else:
            body, shift = '', 0.0
            for t, k in parts:
                target = {'n': 0.0, 'sup': -size * 0.38, 'sub': size * 0.22}[k]
                dy = target - shift
                shift = target
                fsz = f' font-size="{size * 0.7:.2f}"' if k != 'n' else ''
                body += f'<tspan dy="{dy:.2f}"{fsz}>{_esc(t)}</tspan>'
        self.add(f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}"{attrs}{tr}>{body}</text>')

    def render(self, w_mm, h_mm):
        return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w_mm}mm" height="{h_mm}mm" '
                f'viewBox="0 0 {self.w:.2f} {self.h:.2f}" font-family="{_esc(self.font)}">\n'
                f'<rect width="{self.w:.2f}" height="{self.h:.2f}" fill="#FFFFFF"/>\n'
                + '\n'.join(self.el) + '\n</svg>')


# ---------------------------------------------------------------- 主入口
class PlotError(Exception):
    pass


def resolve_colors(names, o, ctrl=None):
    pal = PALETTES.get(o['palette'], PALETTES['nature'])[1]
    out, k = {}, 0
    for g in names:
        if o['colors'].get(g):
            out[g] = o['colors'][g]
            if not (o['ctrl_gray'] and g == ctrl):
                k += 1
            continue
        if o['ctrl_gray'] and g == ctrl:
            out[g] = CTRL_GRAY
            continue
        out[g] = pal[k % len(pal)]
        k += 1
    return out


def make_plot(analysis, opts, gsets=None):
    """analysis: qpcr_analysis.analyze 的结果。opts: 作图参数。返回 SVG 字符串。"""
    o = dict(DEFAULTS)
    o.update({k: v for k, v in (opts or {}).items() if v is not None})
    o['colors'] = o.get('colors') or {}
    if not analysis.get('targets'):
        raise PlotError('还没有可作图的结果。请先在「设置」里选好内参与归一组。')
    groups = [g for g in (o['groups'] or analysis['groups_order']) if g in analysis['groups_order']]
    if not groups:
        groups = analysis['groups_order']
    o['_groups'] = groups
    fn = {'bar': _plot_single, 'box': _plot_single, 'dot': _plot_single,
          'grouped': _plot_grouped, 'heatmap': _plot_heatmap, 'line': _plot_line}.get(o['template'])
    if fn is None:
        raise PlotError(f'未知模板 {o["template"]}')
    return fn(analysis, o, gsets or [])


def _canvas(o):
    W = float(o['width_mm']) * PT_PER_MM
    H = float(o['height_mm']) * PT_PER_MM
    fs = float(o['font_pt'])
    return Svg(W, H, FONTS.get(o['font'], FONTS['arial']), fs), W, H, fs


def _gene_values(gd, g, o):
    """取一个基因某组的 (点, 均值, 误差)，按 y 模式（RQ 或 log2FC）。"""
    d = gd['groups'].get(g)
    if not d:
        return None
    pts = d['log2_points'] if o['y'] == 'log2' else d['rqs']
    if not pts:
        return None
    m = q.mean(pts)
    if len(pts) >= 2:
        e = q.stdev(pts) if o['err'] == 'sd' else q.sem(pts)
    else:
        e = 0.0
    return pts, m, e


def _ylabel(o):
    if o['ylabel']:
        return o['ylabel']
    return 'log_{2} fold change' if o['y'] == 'log2' else 'Relative expression (2^{−ΔΔCt})'


def _annot_for(stats, o, n_groups):
    """返回 (letters, stars)，二者只取其一。"""
    if not stats or stats.get('method') is None or o['annot'] == 'none':
        return None, None
    a = o['annot']
    if a == 'auto':
        a = 'star' if n_groups <= 2 else 'letter'
    if not stats['letters']:  # Dunnett 只比对照，没有字母法
        a = 'star'
    return (stats['letters'], None) if a == 'letter' else (None, stats['stars_vs_control'])


def _y_axis(svg, x0, top, bottom, ticks, y2p, fs, label, o):
    for t in ticks:
        yp = y2p(t)
        svg.line(x0 - 2.2, yp, x0, yp)
        svg.text(x0 - 3.4, yp + fs * 0.35, fmt_tick(t), anchor='end', fill=INK)
    svg.line(x0, y2p(ticks[0]), x0, y2p(ticks[-1]), cap='square')
    tw = max(text_width(fmt_tick(t), fs) for t in ticks)
    lx = x0 - 3.4 - tw - fs * 0.9
    svg.text(lx, (top + bottom) / 2, label, anchor='middle', rotate=-90)


def _left_margin(ticks, fs):
    return max(text_width(fmt_tick(t), fs) for t in ticks) + fs * 2.6 + 6


def _x_labels(svg, centers, labels, base_y, band, fs, italic=False):
    """组名：放得下就水平居中，放不下就 45° 斜排。返回占用高度。"""
    widest = max((text_width(l, fs) for l in labels), default=0)
    if widest <= band * 0.95:
        for cx, l in zip(centers, labels):
            svg.text(cx, base_y + fs * 1.25, l, italic=italic)
        return fs * 1.6
    for cx, l in zip(centers, labels):
        svg.text(cx + fs * 0.3, base_y + fs * 0.9, l, anchor='end', rotate=-45, italic=italic)
    return widest * 0.72 + fs * 1.2


def _x_label_height(labels, band, fs):
    widest = max((text_width(l, fs) for l in labels), default=0)
    return fs * 1.6 if widest <= band * 0.95 else widest * 0.72 + fs * 1.2


def _y_range(vals, o, head_frac):
    lo = min(vals + [0.0])
    hi = max(vals + [0.0])
    if hi == lo:
        hi = lo + 1
    span = hi - lo
    if o['y'] != 'log2':
        lo = 0.0
    hi_need = hi + span * head_frac
    lo_need = lo - (span * 0.04 if lo < 0 else 0)
    ticks = nice_ticks(lo_need, hi_need, 5)
    return ticks


def _fill_style(c, o):
    edge = mix(c, INK, 0.35)
    if o['fill'] == 'light':
        return mix(c, '#FFFFFF', 0.55), edge, 0.6
    if o['fill'] == 'outline':
        return '#FFFFFF', c, 0.9
    return c, edge, 0.4


def _points(svg, pts_xy, c, o, r):
    if o['point'] == 'none':
        return
    for x, y in pts_xy:
        if o['point'] == 'color':
            svg.circle(x, y, r, mix(c, INK, 0.35), '#FFFFFF', 0.3, 0.9)
        elif o['point'] == 'open':
            svg.circle(x, y, r, '#FFFFFF', mix(c, INK, 0.45), 0.5)
        else:
            svg.circle(x, y, r, '#2E2E2E', '#FFFFFF', 0.3, 0.75)


def _errbar(svg, cx, ylo, yhi, cap):
    svg.line(cx, ylo, cx, yhi)
    svg.line(cx - cap, yhi, cx + cap, yhi)
    svg.line(cx - cap, ylo, cx + cap, ylo)


def _signature(svg, W, H, o):
    if o['signature']:
        svg.text(W - 2, H - 2, '自动挡赛车手制作', size=4.5, anchor='end', fill='#BBBBBB')


def _title(svg, W, o, default):
    t = o['title'] if o['title'] else default
    if t:
        svg.text(W / 2, svg.fs * 1.35, t, size=svg.fs * 1.1, weight='600', italic=o['italic'] and not o['title'])
        return svg.fs * 1.9
    return 0


# ---------------------------------------------------------------- 单基因：柱 / 箱 / 点
def _plot_single(a, o, gsets):
    gene = o['gene'] if o['gene'] in a['genes'] else a['targets'][0]
    gd = a['genes'][gene]
    groups = [g for g in o['_groups'] if g in gd['groups']]
    data = {g: _gene_values(gd, g, o) for g in groups}
    groups = [g for g in groups if data[g]]
    if not groups:
        raise PlotError(f'基因 {gene} 没有可用的分组数据。')
    svg, W, H, fs = _canvas(o)
    colors = resolve_colors(groups, o, a['ctrl'])
    letters, stars = _annot_for(gd['stats'], o, len(groups))

    top = 6 + _title(svg, W, o, gene)
    vals = []
    for g in groups:
        pts, m, e = data[g]
        vals += pts + [m + e, m - e]
    if o['template'] == 'box':
        vals = [v for g in groups for v in data[g][0]]
    plot_h_guess = H - top - fs * 3
    head = (fs * 2.0) / max(plot_h_guess, 1) if (letters or stars) else 0.03
    ticks = _y_range(vals, o, head * 1.15)
    left = _left_margin(ticks, fs)
    right = 6
    band = (W - left - right) / len(groups)
    bottom_h = _x_label_height(groups, band, fs) + (fs * 1.6 if o['xlabel'] else 0) + 3
    bottom = H - bottom_h
    y_lo, y_hi = ticks[0], ticks[-1]
    y2p = lambda v: bottom - (v - y_lo) / (y_hi - y_lo) * (bottom - top)
    base = y2p(0 if y_lo <= 0 <= y_hi else y_lo)

    _y_axis(svg, left, top, bottom, ticks, y2p, fs, _ylabel(o), o)
    if y_lo == 0:
        svg.line(left, bottom, W - right, bottom, cap='square')
    else:
        svg.line(left, base, W - right, base, stroke='#9A9A9A', sw=0.4, dash='2,1.5')
    centers = [left + band * (i + 0.5) for i in range(len(groups))]
    bw = min(band * float(o['bar_width']), fs * 6)
    r = max(fs * 0.2, 1.0)

    for g, cx in zip(groups, centers):
        pts, m, e = data[g]
        c = colors[g]
        fill, edge, sw = _fill_style(c, o)
        hi_y = y2p(max(pts))
        if o['template'] == 'bar':
            svg.rect(cx - bw / 2, y2p(m), bw, base - y2p(m), fill, edge, sw)
            if e > 0:
                lo_v = m - e if o['y'] == 'log2' else max(0, m - e)
                _errbar(svg, cx, y2p(lo_v), y2p(m + e), bw * 0.18)
            hi_y = min(hi_y, y2p(m + e), y2p(m))
        elif o['template'] == 'box':
            q1, med, q3 = quantile(pts, .25), quantile(pts, .5), quantile(pts, .75)
            iqr = q3 - q1
            lo_w = min(v for v in pts if v >= q1 - 1.5 * iqr)
            hi_w = max(v for v in pts if v <= q3 + 1.5 * iqr)
            svg.line(cx, y2p(hi_w), cx, y2p(q3), stroke=edge)
            svg.line(cx, y2p(q1), cx, y2p(lo_w), stroke=edge)
            svg.line(cx - bw * 0.2, y2p(hi_w), cx + bw * 0.2, y2p(hi_w), stroke=edge)
            svg.line(cx - bw * 0.2, y2p(lo_w), cx + bw * 0.2, y2p(lo_w), stroke=edge)
            box_fill = mix(c, '#FFFFFF', 0.55) if o['fill'] == 'solid' else fill
            svg.rect(cx - bw / 2, y2p(q3), bw, y2p(q1) - y2p(q3), box_fill, edge, 0.6)
            svg.line(cx - bw / 2, y2p(med), cx + bw / 2, y2p(med), stroke=mix(c, INK, 0.5), sw=1.0)
            hi_y = min(hi_y, y2p(hi_w))
        else:  # dot
            svg.line(cx - bw * 0.36, y2p(m), cx + bw * 0.36, y2p(m), stroke=mix(c, INK, 0.25), sw=1.1)
            if e > 0:
                lo_v = m - e if o['y'] == 'log2' else max(0, m - e)
                _errbar(svg, cx, y2p(lo_v), y2p(m + e), bw * 0.16)
            hi_y = min(hi_y, y2p(m + e))
        pys = [y2p(v) for v in pts]
        span = bw * (0.62 if o['template'] != 'dot' else 0.7)
        offs = jitter(pys, span, r)
        point_o = dict(o)
        if o['template'] == 'dot' and o['point'] == 'dark':
            point_o['point'] = 'color'
        _points(svg, [(cx + dx, py) for dx, py in zip(offs, pys)], c, point_o, r * (1.15 if o['template'] == 'dot' else 1))
        lab = (letters or {}).get(g) or (stars or {}).get(g)
        if lab:
            is_letter = bool(letters)
            svg.text(cx, hi_y - fs * (0.55 if is_letter else 0.35), lab,
                     size=fs * (1.0 if is_letter else 1.15), weight=None if is_letter else '600')

    lab_h = _x_labels(svg, centers, groups, bottom, band, fs)
    if o['xlabel']:
        svg.text((left + W - right) / 2, bottom + lab_h + fs * 1.2, o['xlabel'])
    _signature(svg, W, H, o)
    return svg.render(o['width_mm'], o['height_mm'])


# ---------------------------------------------------------------- 多基因分组柱
def _legend_row(svg, items, colors, x_right, y, fs, kind='box', o=None):
    """右对齐的一行图例，items 为名称列表。返回占用高度。"""
    sw = fs * 0.85
    widths = [sw + 2.5 + text_width(it, fs) + fs * 1.0 for it in items]
    x = x_right - sum(widths)
    for it, w in zip(items, widths):
        c = colors[it]
        if kind == 'line':
            svg.line(x, y - fs * 0.32, x + sw * 1.3, y - fs * 0.32, stroke=c, sw=1.0)
            svg.circle(x + sw * 0.65, y - fs * 0.32, fs * 0.22, c)
            x2 = x + sw * 1.3 + 2.5
        else:
            fill, edge, lw = _fill_style(c, o)
            svg.rect(x, y - sw * 0.85, sw, sw, fill, edge, lw)
            x2 = x + sw + 2.5
        svg.text(x2, y, it, anchor='start')
        x += w + (sw * 0.3 if kind == 'line' else 0)
    return fs * 1.6


def _plot_grouped(a, o, gsets):
    genes = [g for g in (o['genes'] or a['targets']) if g in a['genes']] or a['targets']
    groups = [g for g in o['_groups'] if any(g in a['genes'][t]['groups'] for t in genes)]
    svg, W, H, fs = _canvas(o)
    colors = resolve_colors(groups, o, a['ctrl'])
    top = 6 + _title(svg, W, o, '')
    top += _legend_row(svg, groups, colors, W - 6, top + fs * 0.9, fs, 'box', o)

    vals, cell = [], {}
    for t in genes:
        for g in groups:
            v = _gene_values(a['genes'][t], g, o)
            if v:
                cell[(t, g)] = v
                vals += v[0] + [v[1] + v[2], v[1] - v[2]]
    if not cell:
        raise PlotError('所选基因没有数据。')
    any_annot = o['annot'] != 'none'
    head = (fs * 2.0) / max(H - top - fs * 3, 1) if any_annot else 0.03
    ticks = _y_range(vals, o, head * 1.15)
    left = _left_margin(ticks, fs)
    right = 6
    band = (W - left - right) / len(genes)
    bottom = H - _x_label_height(genes, band, fs) - (fs * 1.6 if o['xlabel'] else 0) - 3
    y_lo, y_hi = ticks[0], ticks[-1]
    y2p = lambda v: bottom - (v - y_lo) / (y_hi - y_lo) * (bottom - top)
    base = y2p(0)
    _y_axis(svg, left, top, bottom, ticks, y2p, fs, _ylabel(o), o)
    if y_lo == 0:
        svg.line(left, bottom, W - right, bottom, cap='square')
    else:
        svg.line(left, base, W - right, base, stroke='#9A9A9A', sw=0.4, dash='2,1.5')

    inner = band * 0.82
    bw = inner / len(groups)
    r = max(fs * 0.17, 0.9)
    centers = []
    for ti, t in enumerate(genes):
        c0 = left + band * (ti + 0.5)
        centers.append(c0)
        gd = a['genes'][t]
        present = [g for g in groups if (t, g) in cell]
        letters, stars = _annot_for(gd['stats'], o, len(present))
        for gi, g in enumerate(groups):
            if (t, g) not in cell:
                continue
            pts, m, e = cell[(t, g)]
            cx = c0 - inner / 2 + bw * (gi + 0.5)
            fill, edge, sw = _fill_style(colors[g], o)
            w = bw * 0.84
            svg.rect(cx - w / 2, y2p(m), w, base - y2p(m), fill, edge, sw)
            if e > 0:
                lo_v = m - e if o['y'] == 'log2' else max(0, m - e)
                _errbar(svg, cx, y2p(lo_v), y2p(m + e), w * 0.2)
            pys = [y2p(v) for v in pts]
            offs = jitter(pys, w * 0.6, r)
            _points(svg, [(cx + dx, py) for dx, py in zip(offs, pys)], colors[g], o, r)
            lab = (letters or {}).get(g) or (stars or {}).get(g)
            if lab:
                hi_y = min(pys + [y2p(m + e), y2p(m)])
                svg.text(cx, hi_y - fs * 0.45, lab, size=fs * (0.9 if letters else 1.05),
                         weight=None if letters else '600')
    lab_h = _x_labels(svg, centers, genes, bottom, band, fs, italic=o['italic'])
    if o['xlabel']:
        svg.text((left + W - right) / 2, bottom + lab_h + fs * 1.2, o['xlabel'])
    _signature(svg, W, H, o)
    return svg.render(o['width_mm'], o['height_mm'])


# ---------------------------------------------------------------- 热图
def _plot_heatmap(a, o, gsets):
    genes = [g for g in (o['genes'] or a['targets']) if g in a['genes']] or a['targets']
    groups = [g for g in o['_groups'] if any(g in a['genes'][t]['groups'] for t in genes)]
    svg, W, H, fs = _canvas(o)
    stops = HEATMAPS.get(o['heat'], HEATMAPS['rdbu'])[1]
    top = 6 + _title(svg, W, o, '')
    vmax = 0.0
    for t in genes:
        for g in groups:
            d = a['genes'][t]['groups'].get(g)
            if d and d['log2fc'] is not None:
                vmax = max(vmax, abs(d['log2fc']))
    # 色条上限：大于 1 取整到整数，否则取到 0.5 的倍数，避免色阶被过度压缩
    vmax = math.ceil(vmax) if vmax > 1 else (math.ceil(vmax * 2) / 2 if vmax > 0 else 1.0)
    left = max(text_width(t, fs) for t in genes) + 8
    cb_w = fs * 0.9
    right = cb_w + text_width(fmt_tick(-vmax), fs) + fs * 3.2
    band = (W - left - right) / len(groups)
    bottom = H - _x_label_height(groups, band, fs) - 3
    ch = (bottom - top) / len(genes)
    for ri, t in enumerate(genes):
        y = top + ch * ri
        svg.text(left - 4, y + ch / 2 + fs * 0.35, t, anchor='end', italic=o['italic'])
        st = a['genes'][t]['stats']
        for ci, g in enumerate(groups):
            x = left + band * ci
            d = a['genes'][t]['groups'].get(g)
            if not d or d['log2fc'] is None:
                svg.rect(x, y, band, ch, '#F2F2F2', '#FFFFFF', 0.8)
                continue
            v = d['log2fc']
            fill = ramp(stops, v / vmax)
            svg.rect(x, y, band, ch, fill, '#FFFFFF', 0.8)
            lum = sum(c * k for c, k in zip(_hex(fill), (0.299, 0.587, 0.114)))
            tc = '#FFFFFF' if lum < 140 else INK
            star = st['stars_vs_control'].get(g, '') if o['annot'] != 'none' else ''
            star = '' if star == 'ns' else star
            if o['heat_values']:
                svg.text(x + band / 2, y + ch / 2 + fs * 0.35, f'{v + 0.0:.2f}'.replace('-0.00', '0.00').replace('-', '−') + star,
                         size=fs * 0.9, fill=tc)
            elif star:
                svg.text(x + band / 2, y + ch / 2 + fs * 0.45, star, size=fs * 1.1, fill=tc, weight='600')
    _x_labels(svg, [left + band * (i + 0.5) for i in range(len(groups))], groups, bottom - 1, band, fs)
    # 色条
    cx0 = W - right + fs * 1.0
    cb_top, cb_bot = max(top + ch * 0.2, top + fs * 2.2), min(bottom, top + max(ch * len(genes) * 0.8, fs * 8))
    n = 40
    for i in range(n):
        t = 1 - 2 * (i + 0.5) / n
        svg.rect(cx0, cb_top + (cb_bot - cb_top) * i / n, cb_w, (cb_bot - cb_top) / n + 0.3, ramp(stops, t))
    for v in (vmax, 0, -vmax):
        yy = cb_top + (cb_bot - cb_top) * (1 - (v + vmax) / (2 * vmax))
        svg.line(cx0 + cb_w, yy, cx0 + cb_w + 1.8, yy)
        svg.text(cx0 + cb_w + 3, yy + fs * 0.33, fmt_tick(v), anchor='start', size=fs * 0.9)
    svg.text(cx0 + cb_w / 2, cb_top - fs * 1.25, 'log_{2}FC', size=fs * 0.9)
    _signature(svg, W, H, o)
    return svg.render(o['width_mm'], o['height_mm'])


# ---------------------------------------------------------------- 时间序列
def _plot_line(a, o, gsets):
    gene = o['gene'] if o['gene'] in a['genes'] else a['targets'][0]
    gd = a['genes'][gene]
    gmap = {s['group']: s for s in gsets}
    groups = [g for g in o['_groups'] if g in gd['groups']]
    series_of = {g: (gmap.get(g, {}).get('series') or '') for g in groups}
    x_of = {g: (gmap.get(g, {}).get('x') or g) for g in groups}
    series = []
    for g in groups:
        if series_of[g] not in series:
            series.append(series_of[g])
    xs_raw = []
    for g in groups:
        if x_of[g] not in xs_raw:
            xs_raw.append(x_of[g])
    numeric = all(_to_float(x) is not None for x in xs_raw) and len(xs_raw) > 1
    if numeric:
        xs_raw.sort(key=_to_float)

    svg, W, H, fs = _canvas(o)
    scolors = resolve_colors([s or gene for s in series], o, None)
    top = 6 + _title(svg, W, o, gene)
    if len(series) > 1:
        top += _legend_row(svg, series, {s: scolors[s] for s in series}, W - 6, top + fs * 0.9, fs, 'line', o)

    vals = []
    data = {}
    for g in groups:
        v = _gene_values(gd, g, o)
        if v:
            data[g] = v
            vals += v[0] + [v[1] + v[2], v[1] - v[2]]
    if not data:
        raise PlotError(f'基因 {gene} 没有可用数据。')
    any_annot = o['annot'] != 'none'
    head = (fs * 2.0) / max(H - top - fs * 3, 1) if any_annot else 0.03
    ticks = _y_range(vals, o, head * 1.15)
    left = _left_margin(ticks, fs) + 2
    right = 8
    if numeric:
        xv = [_to_float(x) for x in xs_raw]
        x_lo, x_hi = min(xv), max(xv)
        pad = (x_hi - x_lo) * 0.06 or 1
        x2p = lambda x: left + (_to_float(x) - x_lo + pad) / (x_hi - x_lo + 2 * pad) * (W - left - right)
        band = (W - left - right) / max(len(xs_raw), 1)
    else:
        band = (W - left - right) / len(xs_raw)
        x2p = lambda x: left + band * (xs_raw.index(x) + 0.5)
    bottom = H - _x_label_height(xs_raw, band, fs) - (fs * 1.6 if o['xlabel'] else 0) - 3
    y_lo, y_hi = ticks[0], ticks[-1]
    y2p = lambda v: bottom - (v - y_lo) / (y_hi - y_lo) * (bottom - top)
    _y_axis(svg, left, top, bottom, ticks, y2p, fs, _ylabel(o), o)
    svg.line(left, bottom, W - right, bottom, cap='square')
    if y_lo < 0:
        svg.line(left, y2p(0), W - right, y2p(0), stroke='#9A9A9A', sw=0.4, dash='2,1.5')
    for x in xs_raw:
        svg.line(x2p(x), bottom, x2p(x), bottom + 2.2)

    r = max(fs * 0.24, 1.1)
    nser = len(series)
    dodge = min(band * 0.12, fs * 0.8) if nser > 1 else 0
    top_at_x = {}
    for si, s in enumerate(series):
        c = scolors[s or gene]
        off = (si - (nser - 1) / 2) * dodge
        pts_line = []
        for x in xs_raw:
            g = next((g for g in groups if series_of[g] == s and x_of[g] == x and g in data), None)
            if not g:
                continue
            pts, m, e = data[g]
            px = x2p(x) + off
            pts_line.append((px, y2p(m)))
            if o['point'] != 'none':
                for v in pts:
                    svg.circle(px, y2p(v), r * 0.7, c, 'none', 0, 0.35)
            if e > 0:
                lo_v = m - e if o['y'] == 'log2' else max(0, m - e)
                _errbar(svg, px, y2p(lo_v), y2p(m + e), r * 1.2)
            top_at_x[x] = min(top_at_x.get(x, 1e9), y2p(m + e), min(y2p(v) for v in pts))
        if len(pts_line) > 1:
            svg.path('M' + ' L'.join(f'{px:.2f},{py:.2f}' for px, py in pts_line), c, 1.0)
        for px, py in pts_line:
            if o['fill'] == 'outline':
                svg.circle(px, py, r * 1.25, '#FFFFFF', c, 0.9)
            else:
                svg.circle(px, py, r * 1.25, c, '#FFFFFF', 0.5)

    # 显著性：多条系列时逐时间点比较各系列（Holm 校正时间点间多重比较），单系列时用整体字母法
    if any_annot:
        if nser == 1:
            st = gd['stats'] if gd['stats']['method'] else {}
            marks = st.get('letters') or {g: s for g, s in (st.get('stars_vs_control') or {}).items() if s}
            for g in groups:
                if g in marks and g in data:
                    svg.text(x2p(x_of[g]), top_at_x[x_of[g]] - fs * 0.6, marks[g], size=fs)
        else:
            sv = a['params']['statval']
            tested = []
            for x in xs_raw:
                sub = {series_of[g]: gd['groups'][g] for g in groups if x_of[g] == x and g in gd['groups']}
                if len(sub) < 2:
                    continue
                res = q.run_stats(sub, value_key=sv)
                if res['method'] is not None:
                    tested.append((x, res['p']))
            for (x, _), padj in zip(tested, q.holm_adjust([p for _, p in tested])):
                lab = q.p_to_stars(padj)
                if lab != 'ns' or o['annot'] in ('star', 'letter'):
                    svg.text(x2p(x), top_at_x[x] - fs * 0.5, lab, size=fs * 1.05, weight='600' if lab != 'ns' else None)

    lab_h = _x_labels(svg, [x2p(x) for x in xs_raw], [str(x) for x in xs_raw], bottom, band, fs)
    if o['xlabel']:
        svg.text((left + W - right) / 2, bottom + lab_h + fs * 1.2, o['xlabel'])
    _signature(svg, W, H, o)
    return svg.render(o['width_mm'], o['height_mm'])


def _to_float(x):
    try:
        return float(str(x).rstrip('dhDHminMIN天小时周 '))
    except ValueError:
        return None


def meta():
    """给前端的模板、配色、预设清单。"""
    return {
        'templates': [{'key': k, 'label': l} for k, l in TEMPLATES],
        'palettes': [{'key': k, 'label': v[0], 'colors': v[1]} for k, v in PALETTES.items()],
        'heatmaps': [{'key': k, 'label': v[0], 'colors': v[1]} for k, v in HEATMAPS.items()],
        'presets': PRESETS, 'defaults': {k: v for k, v in DEFAULTS.items()},
        'fonts': [{'key': k, 'label': {'arial': 'Arial', 'helvetica': 'Helvetica', 'times': 'Times New Roman'}[k]} for k in FONTS],
    }
