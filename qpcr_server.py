# -*- coding: utf-8 -*-
"""
qPCR 数据中心 · 本地服务入口

双击 exe（或 python qpcr_server.py）后在 127.0.0.1 起服务并自动打开浏览器。
只监听本机回环地址，数据不出电脑。浏览器页面关闭 15 分钟后（打包版）自动退出。
"""

import base64
import json
import os
import sys
import tempfile
import threading
import time
import traceback
import urllib.parse
import urllib.request
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import qpcr_core as q
import qpcr_analysis as qa
import qpcr_plots as qp
import qpcr_store as qs

SIGNATURE = 'qpcr-data-center'
VERSION = '2.2'
PORT0 = 8765
IDLE_EXIT = 15 * 60     # 从未连上网页时的兜底退出时间
TAB_STALE = 150         # 标签页超过这么久没心跳视为已关（后台标签的计时器会被浏览器降到每分钟一次）
CLOSE_GRACE = 8         # 最后一个标签页关闭后等这么久再退出，留给刷新页面
FROZEN = getattr(sys, 'frozen', False)
BASE = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
CONFIG = os.path.join(os.path.expanduser('~'), '.qpcr_data_center.json')

state = {'last_ping': time.time(), 'tabs': {}, 'empty_since': None}


def load_config():
    try:
        with open(CONFIG, encoding='utf-8') as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def save_config(cfg):
    with open(CONFIG, 'w', encoding='utf-8') as f:
        json.dump(cfg, f, ensure_ascii=False, indent=1)


def get_store():
    root = load_config().get('root') or qs.default_root()
    os.makedirs(root, exist_ok=True)
    st = state.get('store')
    if st is None or st.root != root:
        st = qs.Store(root)
        state['store'] = st
    return st


def open_path(path):
    if os.name == 'nt':
        os.startfile(path)
    elif sys.platform == 'darwin':
        os.system(f'open "{path}"')
    else:
        os.system(f'xdg-open "{path}" >/dev/null 2>&1 &')


# ---------------------------------------------------------------- API
def api_parse(body):
    """解析上传的 CFX 文件（不入库）。返回孔数据与智能默认。"""
    wells, files = [], []
    tmpdir = tempfile.mkdtemp(prefix='qpcr_')
    try:
        for f in body.get('files', []):
            name = os.path.basename(f['name'])
            path = os.path.join(tmpdir, qs.safe_file(name))
            with open(path, 'wb') as fh:
                fh.write(base64.b64decode(f['data']))
            try:
                ws = q.parse_cfx(path)
                for w in ws:
                    w['file'] = name
                    w['excluded'] = False
                    for k in ('well', 'fluor', 'content'):
                        w[k] = qs._cell_str(w.get(k))
                wells.extend(ws)
                files.append({'name': name, 'n': len(ws),
                              'targets': sorted({w['target'] for w in ws}, key=qa.natural_key),
                              'samples': len({w['sample'] for w in ws if w['sample']}),
                              'error': '' if ws else '未找到有效孔'})
            except Exception as e:
                files.append({'name': name, 'n': 0, 'targets': [], 'samples': 0, 'error': str(e)})
    finally:
        for fn in os.listdir(tmpdir):
            try:
                os.remove(os.path.join(tmpdir, fn))
            except OSError:
                pass
        os.rmdir(tmpdir)
    return {'files': files, 'wells': wells}


def _suggest(wells):
    targets = sorted({w['target'] for w in wells if w['target']}, key=qa.natural_key)
    samples = {w['sample']: q.auto_group_name(w['sample']) for w in wells if w['sample']}
    groups = sorted(set(samples.values()), key=qa.natural_key)
    return {'ref': qa.guess(targets, qa.COMMON_REF) or '',
            'ctrl': qa.guess(groups, qa.COMMON_CTRL) or '',
            'samples': samples}


def api_experiment(params):
    st = get_store()
    p = st.load(params['project'])
    eid = params['id']
    e = p.exp(eid)
    if e is None:
        raise qs.StoreError(f'实验 {eid} 不存在。')
    a = st.analysis(p, eid)
    if p.name in st._disk_dirty:
        st._flush_cache(p.name)
    return {'exp': e, 'wells': p.wells.get(eid, []), 'samples': p.samples.get(eid, {}),
            'gsets': p.gsets.get(eid, []), 'analysis': a}


def _exp_from_body(body):
    e = dict(body.get('exp') or {})
    for k in qs.EXP_KEYS:
        e.setdefault(k, '' if k not in ('plot', 'files') else ({} if k == 'plot' else []))
    e['drop'] = bool(e.get('drop', True))
    return e


def api_analyze(body):
    e = _exp_from_body(body)
    return qa.analyze_cached(e, body.get('wells', []), body.get('samples', {}), body.get('gsets', []))[1]


def api_save(body):
    st = get_store()
    e = _exp_from_body(body)
    files = [(os.path.basename(f['name']), base64.b64decode(f['data'])) for f in body.get('files', [])]
    eid = st.save_experiment(body['project'], e, body.get('wells', []), body.get('samples', {}),
                             body.get('gsets', []), files)
    return {'id': eid}


def api_plot(body):
    st = get_store()
    p = st.load(body['project'])
    eid = body['id']
    a = st.analysis(p, eid)
    svg = qp.make_plot(a, body.get('opts', {}), p.gsets.get(eid, []))
    if body.get('remember'):
        st.update_plot(body['project'], eid, body.get('opts', {}))
    return {'svg': svg}


def api_save_figure(body):
    st = get_store()
    data = body['svg'].encode('utf-8') if body.get('svg') else base64.b64decode(body['png'])
    path = st.save_figure(body['project'], body['filename'], data)
    return {'path': path}


def api_config(body):
    root = (body.get('root') or '').strip()
    cfg = load_config()
    if root:
        os.makedirs(root, exist_ok=True)
        cfg['root'] = root
    else:
        cfg.pop('root', None)
    save_config(cfg)
    return {'root': get_store().root}


def route(method, path, params, body):
    st = get_store
    if path == '/api/ping':
        state['last_ping'] = time.time()
        if params.get('tab'):
            state['tabs'][params['tab']] = time.time()
        return {'ok': True, 'app': SIGNATURE}
    if path == '/api/bye':
        state['tabs'].pop(params.get('tab'), None)
        return {'ok': True}
    if path == '/api/meta':
        return {**qp.meta(), 'root': st().root, 'version': VERSION, 'author': qs.AUTHOR}
    if path == '/api/projects':
        return {'projects': st().list_projects()}
    if path == '/api/experiment':
        return api_experiment(params)
    if method != 'POST':
        raise qs.StoreError('未知请求')
    if path == '/api/project/create':
        p = st().create_project(body.get('name', ''), body.get('desc', ''))
        return {'name': p.name}
    if path == '/api/open_folder':
        target = st().proj_dir(body['project']) if body.get('project') else st().root
        if body.get('sub'):
            target = os.path.join(target, body['sub'])
        os.makedirs(target, exist_ok=True)
        open_path(target)
        return {'ok': True}
    if path == '/api/open_file':
        open_path(st().proj_xlsx(body['project']))
        return {'ok': True}
    if path == '/api/parse':
        r = api_parse(body)
        r['suggest'] = _suggest(r['wells'])
        return r
    if path == '/api/analyze':
        return api_analyze(body)
    if path == '/api/experiment/save':
        return api_save(body)
    if path == '/api/experiment/delete':
        st().delete_experiment(body['project'], body['id'])
        return {'ok': True}
    if path == '/api/plot':
        return api_plot(body)
    if path == '/api/figure/save':
        return api_save_figure(body)
    if path == '/api/config':
        return api_config(body)
    if path == '/api/shutdown':
        threading.Timer(0.3, lambda: os._exit(0)).start()
        return {'ok': True}
    raise qs.StoreError('未知请求')


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def _send(self, code, data, ctype='application/json; charset=utf-8'):
        if not isinstance(data, bytes):
            data = json.dumps(data, ensure_ascii=False).encode('utf-8')
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(data)

    def _handle(self, method):
        u = urllib.parse.urlparse(self.path)
        if not u.path.startswith('/api/'):
            name = 'index.html' if u.path in ('/', '/index.html') else None
            if name is None:
                return self._send(404, b'not found', 'text/plain')
            with open(os.path.join(BASE, 'web', name), 'rb') as f:
                return self._send(200, f.read(), 'text/html; charset=utf-8')
        # 只接受本机页面发起的写请求，防止其它网页跨站调用
        origin = self.headers.get('Origin')
        if method == 'POST' and origin and urllib.parse.urlparse(origin).hostname not in ('127.0.0.1', 'localhost'):
            return self._send(403, {'error': '拒绝跨站请求'})
        params = {k: v[0] for k, v in urllib.parse.parse_qs(u.query).items()}
        body = {}
        if method == 'POST':
            n = int(self.headers.get('Content-Length') or 0)
            if n:
                body = json.loads(self.rfile.read(n).decode('utf-8'))
        state['last_ping'] = time.time()
        try:
            self._send(200, route(method, u.path, params, body))
        except (qs.StoreError, qp.PlotError) as e:
            self._send(400, {'error': str(e)})
        except Exception as e:
            self._send(500, {'error': f'{type(e).__name__}: {e}', 'trace': traceback.format_exc()})

    def do_GET(self):
        self._handle('GET')

    def do_POST(self):
        self._handle('POST')


def _existing_instance(port):
    try:
        with urllib.request.urlopen(f'http://127.0.0.1:{port}/api/ping', timeout=1.5) as r:
            return json.loads(r.read().decode()).get('app') == SIGNATURE
    except Exception:
        return False


def _idle_watch():
    """网页全部关闭后自动退出。从没有网页连上时按 IDLE_EXIT 兜底。"""
    seen_tab = False
    while True:
        time.sleep(2)
        now = time.time()
        tabs = state['tabs']
        for k, t in list(tabs.items()):
            if now - t > TAB_STALE:
                tabs.pop(k, None)
        if tabs:
            seen_tab = True
            state['empty_since'] = None
        elif seen_tab:
            if state['empty_since'] is None:
                state['empty_since'] = now
            elif now - state['empty_since'] >= CLOSE_GRACE:
                os._exit(0)
        if now - state['last_ping'] > IDLE_EXIT:
            os._exit(0)


def main():
    no_browser = '--no-browser' in sys.argv
    srv = None
    for port in range(PORT0, PORT0 + 20):
        if _existing_instance(port):
            if not no_browser:
                webbrowser.open(f'http://127.0.0.1:{port}/')
            return
        try:
            srv = ThreadingHTTPServer(('127.0.0.1', port), Handler)
            break
        except OSError:
            continue
    if srv is None:
        print('找不到可用端口')
        return
    url = f'http://127.0.0.1:{srv.server_address[1]}/'
    print(f'qPCR 数据中心已启动：{url}  （关闭此窗口即退出）')
    if FROZEN:
        threading.Thread(target=_idle_watch, daemon=True).start()
    if not no_browser:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    srv.serve_forever()


if __name__ == '__main__':
    main()
