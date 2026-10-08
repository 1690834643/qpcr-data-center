# -*- coding: utf-8 -*-
"""qPCR 数据中心回归测试。运行：python -m unittest discover -s tests"""

import glob
import os
import shutil
import sys
import tempfile
import unittest
import zipfile
from xml.etree import ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import qpcr_core as q          # noqa: E402
import qpcr_analysis as qa     # noqa: E402
import qpcr_plots as qp        # noqa: E402
import qpcr_store as qs        # noqa: E402

EXAMPLES = os.path.join(os.path.dirname(HERE), '示例数据')
PLATES = sorted(glob.glob(os.path.join(EXAMPLES, '01_*.xlsx')))


def load_plates():
    """示例数据场景一：siRNA 敲低，内参板与目的基因板分开。"""
    wells = []
    for f in PLATES:
        for w in q.parse_cfx(f):
            w['file'] = os.path.basename(f)
            w['excluded'] = False
            wells.append(w)
    return wells


EXP = {'name': 'siRNA 敲低', 'category': '敲低', 'date': '2026-09-12', 'ref': 'GAPDH',
       'ctrl': 'siNC', 'thr': 0.5, 'metric': 'sd', 'drop': True, 'statval': 'dcts',
       'notes': '', 'plot': {}}


class StoreTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.st = qs.Store(self.tmp)
        self.st.create_project('测试项目')
        self.wells = load_plates()
        self.samples = {w['sample']: q.auto_group_name(w['sample']) for w in self.wells if w['sample']}
        self.eid = self.st.save_experiment('测试项目', dict(EXP), self.wells, self.samples, [])

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def fresh(self):
        return qs.Store(self.tmp)

    def test_example_results(self):
        """两块板按样本合并，结果、星号、字母与示例设计一致。"""
        st = self.fresh()
        a = st.analysis(st.load('测试项目'), self.eid)
        expect = {('STAT3', 'siSTAT3'): (0.246, 'b', '***'), ('NFKB1', 'siNFKB1'): (0.133, 'c', '***'),
                  ('NFKB1', 'siSTAT3'): (0.580, 'b', '**'), ('STAT3', 'siNFKB1'): (1.131, 'a', 'ns')}
        for (t, g), (m, letter, star) in expect.items():
            d = a['genes'][t]
            self.assertAlmostEqual(d['groups'][g]['mean'], m, places=3)
            self.assertEqual(d['stats']['letters'][g], letter)
            self.assertEqual(d['stats']['stars_vs_control'][g], star)
        self.assertEqual(a['groups_order'][0], 'siNC')
        self.assertEqual(a['qc_summary']['fixed'], 1)
        self.assertEqual(a['qc_summary']['fail'], 1)
        self.assertNotIn('', {r['sample'] for r in a['qc']})
        self.assertTrue(any('NTC' in x for x in a['warnings']))

    def test_dunnett_option(self):
        """选 Dunnett 后只比各组 vs 对照，不给字母，设置能存回 Excel。"""
        e = dict(EXP, id=self.eid, comp='dunnett')
        self.st.save_experiment('测试项目', e, self.wells, self.samples, [])
        st = self.fresh()
        p = st.load('测试项目')
        self.assertEqual(p.exp(self.eid)['comp'], 'dunnett')
        a = st.analysis(p, self.eid)
        s = a['genes']['STAT3']['stats']
        self.assertEqual(s['method'], 'dunnett')
        self.assertEqual(s['letters'], {})
        self.assertTrue(all(c['group_j'] == 'siNC' for c in s['pairwise']))
        self.assertIn('Dunnett', a['methods_text'])

    def test_excel_resave_roundtrip(self):
        """用户在 Excel 里改了分组与日期后，重新打开按新数据重算。"""
        try:
            import openpyxl
        except ImportError:
            self.skipTest('无 openpyxl')
        path = self.st.proj_xlsx('测试项目')
        wb = openpyxl.load_workbook(path)
        for row in wb['样本分组'].iter_rows(min_row=2):
            if row[1].value == 'siSTAT3-1':
                row[2].value = 'siNC'
        wb['实验目录']['D2'].value = __import__('datetime').date(2026, 9, 13)
        wb.save(path)
        st = self.fresh()
        p = st.load('测试项目')
        self.assertEqual(p.samples[self.eid]['siSTAT3-1'], 'siNC')
        self.assertEqual(p.exp(self.eid)['date'], '2026-09-13')
        self.assertEqual(st.analysis(p, self.eid)['genes']['STAT3']['groups']['siNC']['n'], 5)

    def test_manual_exclusion(self):
        wells = load_plates()
        target = next(i for i, w in enumerate(wells) if w['sample'] == 'siNFKB1-3' and w['target'] == 'NFKB1')
        wells[target]['excluded'] = True
        self.st.save_experiment('测试项目', dict(EXP, id=self.eid), wells, self.samples, [])
        st = self.fresh()
        p = st.load('测试项目')
        self.assertTrue(any(w['excluded'] for w in p.wells[self.eid]))
        r = next(x for x in st.analysis(p, self.eid)['qc'] if x['sample'] == 'siNFKB1-3' and x['target'] == 'NFKB1')
        self.assertEqual(len(r['manual']), 1)
        self.assertEqual(r['n_raw'], 2)

    def test_backup_and_delete(self):
        self.st.delete_experiment('测试项目', self.eid)
        p = self.fresh().load('测试项目')
        self.assertEqual(p.experiments, [])
        self.assertTrue(os.listdir(os.path.join(self.st.proj_dir('测试项目'), '.备份')))

    def test_all_templates(self):
        a = self.st.analysis(self.st.load('测试项目'), self.eid)
        gsets = [{'group': 'siNC', 'x': '0', 'series': 'A'}, {'group': 'siSTAT3', 'x': '24', 'series': 'A'},
                 {'group': 'siNFKB1', 'x': '24', 'series': 'B'}]
        for tpl, _ in qp.TEMPLATES:
            for y in ('rq', 'log2'):
                for preset in qp.PRESETS:
                    o = {'template': tpl, 'y': y, 'palette': preset['palette'], 'fill': preset['fill'],
                         'point': preset['point'], 'ctrl_gray': preset['ctrl_gray'], 'heat': preset['heat']}
                    svg = qp.make_plot(a, o, gsets)
                    ET.fromstring(svg)
                    self.assertIn('width="89mm"', svg)
                    self.assertNotIn('nan', svg.lower())


class ExampleDataTest(unittest.TestCase):
    def test_all_examples_parse(self):
        """六个示例文件（xlsx 与 csv）都能读，三个场景都能算出结果。"""
        cases = [('01_*', 'GAPDH', 'siNC'), ('02_*', 'ACTB', 'Kidney'), ('03_*', 'GAPDH', 'Control_0h')]
        for pat, ref, ctrl in cases:
            files = sorted(glob.glob(os.path.join(EXAMPLES, pat)))
            self.assertEqual(len(files), 2, pat)
            wells = [dict(w, excluded=False) for f in files for w in q.parse_cfx(f)]
            samples = {w['sample']: q.auto_group_name(w['sample']) for w in wells if w['sample']}
            a = qa.analyze({'ref': ref, 'ctrl': ctrl}, wells, samples)
            self.assertTrue(a['targets'], pat)
        self.assertEqual(a['genes']['IL6']['stats']['letters']['LPS_6h'], 'a')


class UnitTest(unittest.TestCase):
    def test_safe_file_keeps_extension(self):
        n = qs.safe_file('E001_' + 'x' * 300 + '.xlsx')
        self.assertTrue(n.endswith('.xlsx'))
        self.assertLessEqual(len(n), 150)
        self.assertEqual(qs.safe_file('a/b:c.csv'), 'a_b_c.csv')

    def test_order_tolerates_text(self):
        o = qa.ordered_groups(['b', 'a', 'CK'], 'CK', [{'group': 'b', 'order': '第一'}, {'group': 'a', 'order': '2'}])
        self.assertEqual(o[0], 'a')

    def test_nice_ticks_negative(self):
        t = qp.nice_ticks(-3.4, 2.1)
        self.assertLessEqual(t[0], -3.4)
        self.assertGreaterEqual(t[-1], 2.1)
        self.assertIn(0, t)

    def test_dunnett_reference(self):
        """Dunnett p 与高精度多元 t 积分（scipy multivariate_t，2e7 次采样）一致。"""
        _, comps = q.dunnett([0.1, -0.4, 0.3, 0.0], [[1.5, 1.1, 1.9], [-0.2, 0.4, 0.1, -0.5, 0.3]])
        self.assertAlmostEqual(comps[0]['t'], 5.547528818117303, places=10)
        self.assertAlmostEqual(comps[0]['p'], 0.00067007, places=6)
        self.assertAlmostEqual(comps[1]['p'], 0.99484692, places=6)

    def test_two_groups_student(self):
        r = q.run_stats({'a': {'dcts': [1.5, 1.1, 1.9]}, 'b': {'dcts': [0.1, -0.4, 0.3, 0.0]}})
        self.assertEqual(r['method'], 't-test')
        self.assertAlmostEqual(r['p'], 0.002203563724022782, places=8)

    def test_holm(self):
        out = q.holm_adjust([0.01, 0.04, 0.03, float('nan'), 0.2])
        for x, y in zip([out[i] for i in (0, 1, 2, 4)], [0.04, 0.09, 0.09, 0.2]):
            self.assertAlmostEqual(x, y)
        self.assertNotEqual(out[3], out[3])

    def test_csv_parse(self):
        ws = q.parse_cfx(sorted(glob.glob(os.path.join(EXAMPLES, '*.csv')))[0])
        self.assertTrue(ws)
        self.assertTrue(all('cq' in w for w in ws))

    _PCRD_XML = '''<?xml version="1.0" encoding="utf-8"?><experimentalData2>
<dyeLayersList><dyeLayer plateName="SYBR" RowsCount="8" ColumnsCount="12" /></dyeLayersList>
<wellSamples>
<wellSample sampleId="" wellSampleType="wcSample" plateIndex="0" geneName="GAPDH" conditionName="CK" />
<wellSample sampleId="" wellSampleType="wcStandard" plateIndex="1" geneName="GAPDH" conditionName="CK" />
<wellSample sampleId="" wellSampleType="wcSample" plateIndex="12" geneName="IL6" conditionName="LPS" />
<wellSample sampleId="" wellSampleType="wcEmpty" plateIndex="13" geneName="" conditionName="" />
</wellSamples>
<dataAnalysisParameters selectedStepNumber="2" selectedWellGroupName="g1" />
<dataAnalysisParam WellGroupGUID="g0" StepNumber="2" wellGroupName="Group 1"><computedWellDataArray>
<computedWellData pIndex="0" thresholdCycle="99.0" isComputed="True" /></computedWellDataArray></dataAnalysisParam>
<dataAnalysisParam WellGroupGUID="g1" StepNumber="4" wellGroupName="All Wells"><computedWellDataArray>
<computedWellData pIndex="0" thresholdCycle="88.0" isComputed="True" /></computedWellDataArray></dataAnalysisParam>
<dataAnalysisParam WellGroupGUID="g1" StepNumber="2" wellGroupName="All Wells"><computedWellDataArray>
<computedWellData pIndex="0" thresholdCycle="21.5" isComputed="True" />
<computedWellData pIndex="1" thresholdCycle="NaN" isComputed="False" />
<computedWellData pIndex="12" thresholdCycle="30.25" isComputed="True" /></computedWellDataArray></dataAnalysisParam>
</experimentalData2>'''

    def test_pcrd_parse(self):
        """合成的 pcrd（ZIP 包 XML 与纯 XML 两种形态）按选中分析读出孔位、基因、样本与 Cq。"""
        tmp = tempfile.mkdtemp()
        try:
            paths = []
            p1 = os.path.join(tmp, 'a.pcrd')
            with zipfile.ZipFile(p1, 'w') as z:
                z.writestr('a.pcrd', self._PCRD_XML)
            paths.append(p1)
            p2 = os.path.join(tmp, 'b.pcrd')
            with open(p2, 'w', encoding='utf-8-sig') as fh:
                fh.write(self._PCRD_XML)
            paths.append(p2)
            for p in paths:
                ws = q.parse_cfx(p)
                by_well = {w['well']: w for w in ws}
                self.assertEqual(sorted(by_well), ['A01', 'A02', 'B01'], p)  # wcEmpty 跳过
                self.assertEqual(by_well['A01']['cq'], 21.5, p)             # 选中分析而非其它组 / 熔解步骤
                self.assertEqual((by_well['A01']['target'], by_well['A01']['sample'], by_well['A01']['content']),
                                 ('GAPDH', 'CK', 'Unkn'), p)
                self.assertEqual(by_well['A01']['fluor'], 'SYBR', p)
                self.assertIsNone(by_well['A02']['cq'], p)                  # NaN / 未计算 -> 缺失
                self.assertEqual(by_well['A02']['content'], 'Std', p)
                self.assertEqual(by_well['B01']['cq'], 30.25, p)            # 行优先编号：12 -> B01
                self.assertEqual((by_well['B01']['target'], by_well['B01']['sample']), ('IL6', 'LPS'), p)
        finally:
            shutil.rmtree(tmp)


if __name__ == '__main__':
    unittest.main()
