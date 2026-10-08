# -*- coding: utf-8 -*-
"""
qPCR 数据中心 · 存储层

每个项目一个文件夹，内含同名 xlsx 总表：
  <数据库根>/<项目>/<项目>.xlsx
  <数据库根>/<项目>/原始文件/   导入时归档的 CFX 原文件
  <数据库根>/<项目>/图/         保存的出图
  <数据库根>/<项目>/.备份/      每次写入前的旧版 xlsx（保留最近 20 份）

xlsx 中「项目信息 / 实验目录 / 原始孔 / 样本分组 / 分组设置」是真值，
其余工作表（QC、ΔΔCt、差异汇总、统计检验）每次保存时重算写出，方便直接在 Excel 里二次使用。
在 Excel 里改了原始孔的 Cq 或分组，保存后重新打开项目即按新数据重算。
"""

import datetime
import json
import os
import re
import shutil
import threading

import qpcr_core as q
import qpcr_analysis as qa

APP_NAME = 'qPCR 数据中心'
AUTHOR = '自动挡赛车手制作'

EXP_COLS = ['实验ID', '实验名称', '实验内容', '日期', '内参基因', '归一组', 'QC阈值',
            '离散指标', '剔除离群', '统计值', '原始文件', '备注', '创建时间', '更新时间', '作图设置']
EXP_KEYS = ['id', 'name', 'category', 'date', 'ref', 'ctrl', 'thr',
            'metric', 'drop', 'statval', 'files', 'notes', 'created', 'updated', 'plot']
WELL_COLS = ['实验ID', '来源文件', '孔位', '荧光', '基因', '类型', '样本', 'Cq', '手动剔除']
WELL_KEYS = ['exp', 'file', 'well', 'fluor', 'target', 'content', 'sample', 'cq', 'excluded']
GSET_COLS = ['实验ID', '分组', '顺序', '系列', '时间点', '颜色']
GSET_KEYS = ['exp', 'group', 'order', 'series', 'x', 'color']

_lock = threading.RLock()


class StoreError(Exception):
    """给用户看的存储错误（中文说明）。"""


def now():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M')


def safe_name(name):
    """项目名转文件夹名：去掉 Windows 非法字符与首尾空白点号。"""
    s = re.sub(r'[\\/:*?"<>|\x00-\x1f]', '_', str(name)).strip().strip('.')
    return s[:60]


def safe_file(name, limit=150):
    """文件名去非法字符，过长时截断主名但保留扩展名。"""
    s = re.sub(r'[\\/:*?"<>|\x00-\x1f]', '_', str(name)).strip()
    stem, dot, ext = s.rpartition('.')
    if not dot:
        stem, ext = s, ''
    keep = limit - len(ext) - 1
    return (stem[:keep].rstrip(' .') + ('.' + ext if ext else '')) or 'file'


def default_root():
    """默认数据库位置：我的文档/qPCR数据库（Windows 取真实「文档」路径，兼容 OneDrive 重定向）。"""
    env = os.environ.get('QPCR_DATA_ROOT')
    if env:
        return env
    docs = None
    if os.name == 'nt':
        try:
            import ctypes
            from ctypes import wintypes
            buf = ctypes.create_unicode_buffer(wintypes.MAX_PATH)
            # CSIDL_PERSONAL = 5（我的文档）
            if ctypes.windll.shell32.SHGetFolderPathW(None, 5, None, 0, buf) == 0:
                docs = buf.value
        except Exception:
            docs = None
    if not docs:
        docs = os.path.join(os.path.expanduser('~'), 'Documents')
    return os.path.join(docs, 'qPCR数据库')


def _cell_str(v):
    if v is None:
        return ''
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v).strip()


def _excel_date(v):
    """Excel 把日期改成序列号时还原为 YYYY-MM-DD。"""
    if isinstance(v, float) and 20000 < v < 80000:
        d = datetime.date(1899, 12, 30) + datetime.timedelta(days=int(v))
        return d.isoformat()
    return _cell_str(v)


def _table(rows, keys, cols):
    """按表头名取列，返回 list[dict]，表头列顺序被用户调整过也能读。"""
    if not rows:
        return []
    header = [_cell_str(x) for x in rows[0]]
    idx = {c: header.index(c) for c in cols if c in header}
    out = []
    for r in rows[1:]:
        if not r or all(x is None or x == '' for x in r):
            continue
        d = {}
        for k, c in zip(keys, cols):
            i = idx.get(c)
            d[k] = r[i] if i is not None and i < len(r) else None
        out.append(d)
    return out


class Project:
    def __init__(self, name, desc='', created=None):
        self.name = name
        self.desc = desc
        self.created = created or now()
        self.updated = self.created
        self.experiments = []      # list[dict]，键见 EXP_KEYS
        self.wells = {}            # exp_id -> list[dict]
        self.samples = {}          # exp_id -> {sample: group}
        self.gsets = {}            # exp_id -> list[dict]

    def exp(self, exp_id):
        return next((e for e in self.experiments if e['id'] == exp_id), None)

    def next_id(self):
        n = 0
        for e in self.experiments:
            m = re.match(r'E(\d+)$', e['id'])
            if m:
                n = max(n, int(m.group(1)))
        return f'E{n + 1:03d}'

    def categories(self):
        return sorted({e['category'] for e in self.experiments if e.get('category')}, key=qa.natural_key)


class Store:
    def __init__(self, root):
        self.root = root
        self._cache = {}    # 项目名 -> (mtime, Project)
        self._disk = {}     # 项目名 -> {输入指纹: 分析结果}
        self._disk_dirty = set()
        self._used_keys = {}

    # ---------- 路径 ----------
    def proj_dir(self, name):
        return os.path.join(self.root, safe_name(name))

    def proj_xlsx(self, name):
        return os.path.join(self.proj_dir(name), safe_name(name) + '.xlsx')

    # ---------- 项目列表 ----------
    def list_projects(self):
        out = []
        if not os.path.isdir(self.root):
            return out
        for d in sorted(os.listdir(self.root), key=qa.natural_key):
            x = os.path.join(self.root, d, d + '.xlsx')
            if not os.path.isfile(x):
                continue
            try:
                p = self.load(d)
                out.append({'name': p.name, 'desc': p.desc, 'n_exp': len(p.experiments),
                            'updated': p.updated, 'categories': p.categories(),
                            'experiments': [self._exp_card(p, e) for e in p.experiments]})
            except Exception as e:
                out.append({'name': d, 'desc': '', 'n_exp': 0, 'updated': '',
                            'categories': [], 'experiments': [], 'error': f'读取失败：{e}'})
        return out

    @staticmethod
    def _exp_card(p, e):
        ws = p.wells.get(e['id'], [])
        groups = set(p.samples.get(e['id'], {}).values())
        return {'id': e['id'], 'name': e['name'], 'category': e.get('category', ''),
                'date': e.get('date', ''), 'ref': e.get('ref', ''), 'notes': e.get('notes', ''),
                'targets': sorted({w['target'] for w in ws if w['target']}, key=qa.natural_key),
                'n_samples': len({w['sample'] for w in ws if w['sample']}),
                'n_groups': len(groups), 'n_wells': len(ws), 'updated': e.get('updated', '')}

    def create_project(self, name, desc=''):
        name = safe_name(name)
        if not name:
            raise StoreError('项目名不能为空。')
        with _lock:
            if os.path.exists(self.proj_xlsx(name)):
                raise StoreError(f'项目「{name}」已存在。')
            os.makedirs(os.path.join(self.proj_dir(name), '原始文件'), exist_ok=True)
            os.makedirs(os.path.join(self.proj_dir(name), '图'), exist_ok=True)
            p = Project(name, desc)
            self.save(p)
            return p

    # ---------- 读 ----------
    def load(self, name):
        path = self.proj_xlsx(name)
        if not os.path.isfile(path):
            raise StoreError(f'找不到项目「{name}」。')
        mt = os.path.getmtime(path)
        hit = self._cache.get(name)
        if hit and hit[0] == mt:
            return hit[1]
        sheets = q.read_xlsx_sheets(path)
        info = {}
        for r in sheets.get('项目信息', [])[1:]:
            if r and len(r) >= 2:
                info[_cell_str(r[0])] = _cell_str(r[1])
        p = Project(info.get('项目名称') or name, info.get('项目描述', ''), info.get('创建时间'))
        p.updated = info.get('更新时间', p.created)
        for e in _table(sheets.get('实验目录', []), EXP_KEYS, EXP_COLS):
            e = {k: (_excel_date(v) if k == 'date' else _cell_str(v)) for k, v in e.items()}
            if not e['id']:
                continue
            try:
                e['plot'] = json.loads(e['plot']) if e['plot'] else {}
            except ValueError:
                e['plot'] = {}
            e['files'] = [f for f in e['files'].split(' | ') if f] if e['files'] else []
            e['drop'] = e['drop'] in ('是', 'True', 'true', '1')
            e['metric'] = 'range' if e['metric'] in ('极差', 'range') else 'sd'
            e['statval'] = 'rqs' if e['statval'] in ('RQ', 'rqs') else 'dcts'
            p.experiments.append(e)
        for w in _table(sheets.get('原始孔', []), WELL_KEYS, WELL_COLS):
            eid = _cell_str(w['exp'])
            cq = w['cq']
            if isinstance(cq, str):
                try:
                    cq = float(cq)
                except ValueError:
                    cq = None
            if isinstance(cq, float) and cq != cq:
                cq = None
            p.wells.setdefault(eid, []).append({
                'file': _cell_str(w['file']), 'well': _cell_str(w['well']),
                'fluor': _cell_str(w['fluor']), 'target': _cell_str(w['target']),
                'content': _cell_str(w['content']), 'sample': _cell_str(w['sample']),
                'cq': cq, 'excluded': _cell_str(w['excluded']) in ('是', '1', 'True', 'true'),
            })
        for r in _table(sheets.get('样本分组', []), ['exp', 'sample', 'group'], ['实验ID', '样本', '分组']):
            eid = _cell_str(r['exp'])
            s = _cell_str(r['sample'])
            if eid and s:
                p.samples.setdefault(eid, {})[s] = _cell_str(r['group']) or q.auto_group_name(s)
        for r in _table(sheets.get('分组设置', []), GSET_KEYS, GSET_COLS):
            eid = _cell_str(r['exp'])
            g = _cell_str(r['group'])
            if eid and g:
                o = r['order']
                p.gsets.setdefault(eid, []).append({
                    'group': g, 'order': o if isinstance(o, float) else (_cell_str(o) or None),
                    'series': _cell_str(r['series']), 'x': _cell_str(r['x']),
                    'color': _cell_str(r['color'])})
        self._cache[name] = (mt, p)
        return p

    def analysis(self, p, eid):
        e = p.exp(eid)
        if e is None:
            raise StoreError(f'项目「{p.name}」里没有实验 {eid}。')
        wells = p.wells.get(eid, [])
        samples = p.samples.get(eid, {})
        for w in wells:
            if w['sample'] and w['sample'] not in samples:
                samples[w['sample']] = q.auto_group_name(w['sample'])
        disk = self._disk_cache(p.name)
        key, res = qa.analyze_cached(e, wells, samples, p.gsets.get(eid), disk)
        if key not in disk:
            disk[key] = res
            self._disk_dirty.add(p.name)
        self._used_keys.setdefault(p.name, set()).add(key)
        return res

    def _cache_path(self, name):
        return os.path.join(self.proj_dir(name), '.分析缓存.json')

    def _disk_cache(self, name):
        if name not in self._disk:
            try:
                with open(self._cache_path(name), encoding='utf-8') as f:
                    self._disk[name] = json.load(f)
            except (OSError, ValueError):
                self._disk[name] = {}
        return self._disk[name]

    def _flush_cache(self, name, keep=None):
        """写回项目分析缓存，只保留当前各实验用到的条目。"""
        d = self._disk_cache(name)
        if keep is not None:
            for k in [k for k in d if k not in keep]:
                del d[k]
        try:
            with open(self._cache_path(name), 'w', encoding='utf-8') as f:
                json.dump(d, f, ensure_ascii=False)
        except OSError:
            pass
        self._disk_dirty.discard(name)

    # ---------- 写 ----------
    def save_experiment(self, name, exp, wells, samples, gsets, files=None):
        """新建或更新实验。exp 无 id 时新建。files: list[(文件名, bytes)] 原始文件归档。"""
        with _lock:
            p = self.load(name)
            if not exp.get('id'):
                exp['id'] = p.next_id()
                exp['created'] = now()
                p.experiments.append(exp)
            else:
                old = p.exp(exp['id'])
                if old is None:
                    raise StoreError(f'实验 {exp["id"]} 不存在。')
                exp['created'] = old.get('created', now())
                if not exp.get('files'):
                    exp['files'] = old.get('files', [])
                p.experiments[p.experiments.index(old)] = exp
            exp['updated'] = now()
            if files:
                raw_dir = os.path.join(self.proj_dir(name), '原始文件')
                os.makedirs(raw_dir, exist_ok=True)
                saved = []
                for fname, data in files:
                    fn = safe_file(f'{exp["id"]}_{fname}')
                    with open(os.path.join(raw_dir, fn), 'wb') as f:
                        f.write(data)
                    saved.append(fname)
                exp['files'] = [f for f in (exp.get('files') or []) if f not in saved] + saved
            p.wells[exp['id']] = wells
            p.samples[exp['id']] = samples
            p.gsets[exp['id']] = gsets
            self.save(p)
            return exp['id']

    def update_plot(self, name, eid, plot):
        with _lock:
            p = self.load(name)
            e = p.exp(eid)
            if e is None:
                raise StoreError(f'实验 {eid} 不存在。')
            e['plot'] = plot
            self.save(p)

    def delete_experiment(self, name, eid):
        with _lock:
            p = self.load(name)
            p.experiments = [e for e in p.experiments if e['id'] != eid]
            for d in (p.wells, p.samples, p.gsets):
                d.pop(eid, None)
            self.save(p)

    def save_figure(self, name, filename, data):
        d = os.path.join(self.proj_dir(name), '图')
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, safe_file(filename))
        with open(path, 'wb') as f:
            f.write(data)
        return path

    def save(self, p):
        path = self.proj_xlsx(p.name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        p.updated = now()
        self._used_keys[p.name] = set()
        sheets = self._build_sheets(p)
        self._flush_cache(p.name, keep=self._used_keys[p.name])
        if os.path.isfile(path):
            self._backup(p.name, path)
        tmp = path + '.tmp'
        q.write_xlsx(tmp, sheets)
        try:
            os.replace(tmp, path)
        except PermissionError:
            os.remove(tmp)
            raise StoreError(f'无法写入「{os.path.basename(path)}」。它可能正在 Excel / WPS 中打开，请先关闭再保存。')
        self._cache[p.name] = (os.path.getmtime(path), p)

    def _backup(self, name, path):
        bdir = os.path.join(self.proj_dir(name), '.备份')
        os.makedirs(bdir, exist_ok=True)
        stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        shutil.copy2(path, os.path.join(bdir, f'{safe_name(name)}_{stamp}.xlsx'))
        olds = sorted(f for f in os.listdir(bdir) if f.endswith('.xlsx'))
        for f in olds[:-20]:
            try:
                os.remove(os.path.join(bdir, f))
            except OSError:
                pass

    # ---------- 组装 xlsx ----------
    def _build_sheets(self, p):
        H = lambda cols: [(c, q.S_HEADER) for c in cols]
        opt = lambda widths: {'freeze_header': True, 'autofilter': True, 'col_widths': widths}

        info = [H(['项目', '内容']),
                ['项目名称', p.name], ['项目描述', p.desc], ['创建时间', p.created],
                ['更新时间', p.updated], ['实验数', len(p.experiments)],
                ['生成程序', f'{APP_NAME}（{AUTHOR}）']]

        exp_rows = [H(EXP_COLS)]
        well_rows = [H(WELL_COLS)]
        samp_rows = [H(['实验ID', '样本', '分组'])]
        gset_rows = [H(GSET_COLS)]
        qc_rows = [H(['实验ID', '实验名称', '样本', '分组', '基因', '孔数', 'Cq均值', 'SD', '极差',
                      '采用Cq', '自动剔除', '手动剔除', 'QC判定'])]
        res_rows = [H(['实验ID', '实验名称', '样本', '分组', '基因', '内参Cq', '目的Cq',
                       'ΔCt', 'ΔΔCt', 'RQ(2^-ΔΔCt)', 'log2FC'])]
        sum_rows = [H(['实验ID', '实验名称', '实验内容', '基因', '分组', 'n', 'RQ均值', 'SEM', 'SD',
                       'log2FC', 'p(vs对照)', '显著性', '字母', '检验方法'])]
        stat_rows = [H(['实验ID', '实验名称', '基因', '方法', '总体p', '比较', 'p', '显著性'])]
        method_name = {'t-test': 'Welch t 检验', 'anova': '单因素 ANOVA + Tukey HSD'}

        for e in p.experiments:
            eid = e['id']
            exp_rows.append([
                eid, e.get('name', ''), e.get('category', ''), e.get('date', ''),
                e.get('ref', ''), e.get('ctrl', ''), _num_or_str(e.get('thr') or 0.5),
                '极差' if e.get('metric') == 'range' else 'SD',
                '是' if e.get('drop', True) else '否',
                'RQ' if e.get('statval') == 'rqs' else 'ΔCt',
                ' | '.join(e.get('files') or []), e.get('notes', ''),
                e.get('created', ''), e.get('updated', ''),
                json.dumps(e.get('plot') or {}, ensure_ascii=False)])
            for w in p.wells.get(eid, []):
                well_rows.append([eid, w.get('file', ''), w.get('well', ''), w.get('fluor', ''),
                                  w.get('target', ''), w.get('content', ''), w.get('sample', ''),
                                  (round(w['cq'], 4), q.S_NUM2) if w.get('cq') is not None else '',
                                  ('是', q.S_WARN) if w.get('excluded') else ''])
            for s, g in sorted(p.samples.get(eid, {}).items(), key=lambda kv: qa.natural_key(kv[0])):
                samp_rows.append([eid, s, g])
            for g in p.gsets.get(eid, []):
                gset_rows.append([eid, g['group'], _num_or_str(g.get('order')), g.get('series', ''),
                                  g.get('x', ''), g.get('color', '')])

            try:
                a = self.analysis(p, eid)
            except Exception:
                continue
            ename = e.get('name', '')
            st_style = {'pass': q.S_PASS, 'single': q.S_PASS, 'fixed': q.S_WARN, 'fail': q.S_FAIL}
            for r in a['qc']:
                qc_rows.append([eid, ename, r['sample'], r['group'], r['target'], r['n_raw'],
                                (round(r['mean'], 3), q.S_NUM2) if r['mean'] is not None else '',
                                (round(r['sd'] or 0, 4), q.S_NUM4), (round(r['range'] or 0, 4), q.S_NUM4),
                                '/'.join(f'{c:.3f}' for c in r['used']),
                                '/'.join(f'{c:.3f}' for c in r['outliers']),
                                '/'.join(f'{c:.3f}' for c in r['manual']),
                                (r['reason'], st_style[r['status']])])
            for v in a['per_sample']:
                rq = v['rq']
                res_rows.append([eid, ename, v['sample'], v['group'], v['target'],
                                 _r(v['ref_cq'], 3, q.S_NUM2), _r(v['tgt_cq'], 3, q.S_NUM2),
                                 _r(v['dct'], 4, q.S_NUM4), _r(v['ddct'], 4, q.S_NUM4),
                                 _r(rq, 4, q.S_NUM4),
                                 _r(-v['ddct'], 4, q.S_NUM4) if v['ddct'] is not None else ''])
            for t in a['targets']:
                gd = a['genes'][t]
                st = gd['stats']
                for g in a['groups_order']:
                    d = gd['groups'].get(g)
                    if not d:
                        continue
                    sum_rows.append([eid, ename, e.get('category', ''), t, g, d['n'],
                                     _r(d['mean'], 4, q.S_NUM4), _r(d['sem'], 4, q.S_NUM4),
                                     _r(d['sd'], 4, q.S_NUM4), _r(d['log2fc'], 4, q.S_NUM4),
                                     _r(st['p_vs_control'].get(g), 5, q.S_NUM4) if g != a['ctrl'] else '对照',
                                     st['stars_vs_control'].get(g, ''), st['letters'].get(g, ''),
                                     method_name.get(st['method'], '未检验')])
                first = True
                for c in st['pairwise'] or [None]:
                    stat_rows.append([
                        eid if first else '', ename if first else '', t if first else '',
                        method_name.get(st['method'], '未检验') if first else '',
                        _r(st['p'], 5, q.S_NUM4) if first else '',
                        f'{c["group_i"]} vs {c["group_j"]}' if c else '',
                        _r(c['p'], 5, q.S_NUM4) if c else '', c['stars'] if c else ''])
                    first = False

        notes = [H(['说明', '内容']),
                 ['真值表', '项目信息 / 实验目录 / 原始孔 / 样本分组 / 分组设置。修改这些表后重新打开项目即按新数据重算。'],
                 ['派生表', '技术重复QC / ΔΔCt结果 / 差异汇总 / 统计检验。每次保存自动重算，手改会被覆盖。'],
                 ['ΔΔCt', 'ΔCt = Cq(目的) − Cq(内参)，ΔΔCt = ΔCt − mean(ΔCt 归一组)，RQ = 2^(−ΔΔCt)'],
                 ['log2FC', 'log2FC = −mean(ΔΔCt)，即组内 RQ 几何均值的 log2'],
                 ['技术重复QC', '同一样本×基因的复孔 SD（或极差）超过阈值判为离散，可开启三选二自动剔除离群孔'],
                 ['显著性', '2 组用 Welch t 检验，≥3 组用单因素 ANOVA + Tukey HSD，字母法中 a 给 RQ 最高组'],
                 ['星号', '* p<0.05，** p<0.01，*** p<0.001，ns 不显著'],
                 ['生成程序', f'{APP_NAME}（{AUTHOR}）']]

        return [
            ('项目信息', info, opt([14, 60])),
            ('实验目录', exp_rows, opt([8, 22, 14, 12, 10, 12, 8, 8, 8, 8, 40, 30, 16, 16, 30])),
            ('原始孔', well_rows, opt([8, 36, 7, 7, 12, 7, 16, 9, 9])),
            ('样本分组', samp_rows, opt([8, 18, 16])),
            ('分组设置', gset_rows, opt([8, 16, 7, 12, 10, 10])),
            ('技术重复QC', qc_rows, opt([8, 18, 16, 14, 12, 6, 9, 8, 8, 24, 10, 10, 16])),
            ('ΔΔCt结果', res_rows, opt([8, 18, 16, 14, 12, 9, 9, 9, 9, 12, 9])),
            ('差异汇总', sum_rows, opt([8, 18, 12, 12, 14, 5, 10, 9, 9, 9, 11, 8, 6, 22])),
            ('统计检验', stat_rows, opt([8, 18, 12, 22, 10, 26, 10, 8])),
            ('说明', notes, {'col_widths': [14, 90]}),
        ]


def _num_or_str(v):
    if v is None or v == '':
        return ''
    try:
        f = float(v)
        return int(f) if f.is_integer() else f
    except (TypeError, ValueError):
        return str(v)


def _r(v, nd, style):
    if v is None:
        return ''
    return (round(v, nd), style)
