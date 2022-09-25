""" Code for running the C++ versino of the planner and testing it
"""
from trajectory import Joint, Block, Segment, SegmentList
from collections import deque
import subprocess
import os
from pathlib import Path

precision_map = {
    'v': 0,
    'x': 0,
    't': 3,
    'd': 0,
}

def report_planner_diffs(sl1, sl2):
    for seg_n, seg_diffs in compare_planner(sl1, sl2):
        print("=====", seg_n)
        for block_n, block_diffs in seg_diffs:
            print("Block: ", block_n)
            for var, b1, b2 in block_diffs:
                print("   ", var, abs(b1-b2), b1, b2)

def compare_blocks(a, b):
    diffs = []
    for k, v in a.items():
        if k not in ('_tag', '_type', 'replans', 'reductions', 'joint', 'segment', 'memo'):
            precision = precision_map.get(k[0], 0)
            if v is not None and b[k] is not None and round(v, precision) != round(b[k], precision):
                diffs.append((k, v, b[k]))

    return diffs


def compare_planner(sl1, sl2):
    diffs = []
    for i, (s1, s2) in enumerate(zip(sl1.segments, sl2.segments)):
        seg_diffs = compare_seg(s1, s2)
        if seg_diffs:
            diffs.append((i, seg_diffs))

    return diffs

def compare_seg(s1, s2):
    sdiffs = []
    for i, b1 in enumerate(s1.blocks):
        b2 = s2.blocks[i]
        diffs = compare_blocks(b1.asdict(), b2.asdict())
        if diffs:
            sdiffs.append((i, diffs))

    return sdiffs


class TestPlanner:
    """Run the unit test program, looking for tests that output JSON"""
    def __init__(self, test_dir):
        from pathlib import Path
        self.test_dir = Path(test_dir)
        self.exe_path = self.test_dir.joinpath('test_planner')

    def make(self):

        try:
            pwd = os.getcwd()
            os.chdir(self.test_dir)
            os.system('make')

        finally:
            os.chdir(pwd)

    def run(self):

        import subprocess
        import json
        o = {}

        r = subprocess.check_output([self.exe_path, r'[json]'], timeout=1)
        r = r.decode('utf8')

        for l in r.splitlines():
            if l.startswith('JSON'):
                d = json.loads(l[4:])
                o[d['test']] = d

        return o

    def load_segment(self, test_name):

        d = self.run()

        j = d[test_name]['output']

        blocks = []
        joints = []
        for jd in j['joints']:
            jd = dict(**jd)
            del jd['_type']
            joints.append(Joint(**jd))

        for joint, bd in zip(joints, j['blocks']):
            bd = dict(**bd)
            del bd['_type']

            blocks.append(Block(**bd, joint=joint))

        s = Segment(0, joints, blocks)
        s.move =j['move']

        return s

class CPPPlanner:
    """Run the program that recieves std input describing moves, runs it through the planner and
    return JSON results"""

    def __init__(self, test_dir):
        self.test_dir = Path(test_dir)
        self.exe_path = self.test_dir.joinpath('planner')

    def make(self):

        try:
            pwd = os.getcwd()
            os.chdir(self.test_dir)
            os.system('make')

        finally:
            os.chdir(pwd)

    def make_input(self, joints, moves):
        inpt = f"{len(joints)}\n"

        for j in joints:
            inpt += f"{j.v_max} {j.a_max}\n"

        for m in moves:
            inpt += ' '.join(str(e) for e in m) + '\n'

        return inpt


    def _run_planner(self, joints, moves, options):
        import json

        inpt = self.make_input(joints, moves)

        with open(self.test_dir.joinpath('planner_input.txt'), 'w') as f:
            f.write(inpt)

        proc = subprocess.Popen([str(self.exe_path),options], text=True, encoding='utf-8',
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        try:
            outs, errs = proc.communicate(input=inpt, timeout=2)
        except subprocess.TimeoutExpired:
            print("Timedout")
            proc.kill()
            outs, errs = proc.communicate()

        if 'j' in options: # JSON output
            for line in outs.splitlines():
                try:
                    return json.loads(line)
                except json.JSONDecodeError as e:
                    print("| ",line)
        elif 's' in options: # stepper output
            steps = []
            for line in outs.splitlines():
                t, *d = line.split()
                steps.append([float(t)]+list(map(int, d)))
            return steps


    def run_planner(self, joints, moves):
        return self._run_planner(joints, moves, '-pj')

    def run_stepper(self, joints, moves):
        return self._run_planner(joints, moves, '-s')

    def planner(self, joints, moves):

        d = self.run_planner(joints, moves)

        sl = SegmentList(joints)

        joints = [Joint(n=i, v_max=j['v_max'],a_max=j['a_max']) for i,j in enumerate(d['joints'])]
        for seg_n, sd in  enumerate(d['segments']):
            blocks = []
            for joint, bd in zip(joints, sd['blocks']):
                del bd['_type']
                b = Block(**bd, joint=joint)
                blocks.append(b)
            sl.segments.append(Segment(seg_n, joints, blocks))

        if "_time" in d:
            print("CPP Time: ", round(d["_time"]), "μs",
                  round(d["_time"]/(len(sl.segments)*len(joints))), "μs per block")

        return sl

    def compare_planner(self, joints, moves, report=True):
        import time

        sl_p = SegmentList(joints)
        start = time.time()
        for move in moves:
            sl_p.move(move)
        t_diff = time.time()-start
        print("Pyp Time: ", round(t_diff*1_000_000) , "μs",
              round(t_diff*1_000_000/(len(sl_p.segments)*len(joints))), "μs per block")

        self.make()
        sl_c = self.planner(joints, moves)

        if report:
            report_planner_diffs(sl_p, sl_c)

        return sl_p, sl_c


class TestStepper:
    """Run the test programm that generates steps"""

    def __init__(self, test_dir):
        self.test_dir = Path(test_dir)
        self.exe_path = self.test_dir.joinpath('stepper')

    def make(self):

        try:
            pwd = os.getcwd()
            os.chdir(self.test_dir)
            os.system('make')

        finally:
            os.chdir(pwd)

    def make_input(self, sl):
        from itertools import chain
        inpt = ""
        for s in sl.segments:
            sb = [b.stepper_blocks() for b in s]
            for e in zip(*sb):
                inpt += " ".join(str(i) for i in list(chain(*e)))+'\n'

        return inpt


    def run_stepper(self, sl):
        import json

        self.make()

        inpt = self.make_input(sl)

        with open(self.test_dir.joinpath('stepper_input.txt'), 'w') as f:
            f.write(inpt)

        proc = subprocess.Popen(str(self.exe_path), text=True, encoding='utf-8',
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        try:
            outs, errs = proc.communicate(input=inpt, timeout=2)
        except subprocess.TimeoutExpired:
            print("Timedout")
            proc.kill()
            outs, errs = proc.communicate()

        for line in outs.splitlines():
            print(line)

