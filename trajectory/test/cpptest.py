""" Code for running the C++ versino of the planner and testing it
"""
from trajectory import Joint, ACDBlock, Segment, SegmentList
from collections import deque
import subprocess
import os
from pathlib import Path

class CPPPlanner:

    def __init__(self, test_dir):
        self.test_dir = Path(test_dir)
        self.exe_path = self.test_dir.joinpath('trj_planner')

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


    def run_planner(self, joints, moves):
        import json

        inpt = self.make_input(joints, moves)

        proc = subprocess.Popen(str(self.exe_path), text=True, encoding='utf-8',
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        try:
            outs, errs = proc.communicate(input=inpt, timeout=2)
        except subprocess.TimeoutExpired:
            print("Timedout")
            proc.kill()
            outs, errs = proc.communicate()

        try:
            return json.loads(outs)
        except json.JSONDecodeError as e:
            print("ERROR", e)
            print(outs)

    def planner(self, joints, moves):

        d = self.run_planner(joints, moves)

        sl = SegmentList(joints)

        joints = [Joint(n=i, v_max=j['v_max'],a_max=j['a_max']) for i,j in enumerate(d['joints'])]
        for seg_n, sd in  enumerate(d['segments']):
            blocks = []
            for bd in sd['blocks']:
                del bd['_type']
                b = ACDBlock(**bd)
                blocks.append(b)
            sl.segments.append(Segment(seg_n, blocks, joints))

        return sl

