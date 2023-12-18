# -*- coding:utf-8 -*-
import os
import argparse
from . import config
from .app import main

parser = argparse.ArgumentParser(description='Sān-D Molecule')
parser.add_argument('files', type=str, nargs='*',
                    help='open files by SDM')

args = parser.parse_args()

cwd = os.getcwd()
try:
    os.chdir(config['PATH']['module'])
    if len(args.files) > 0:
        main(args.files)
    else:
        main()
finally:
    os.chdir(cwd)
