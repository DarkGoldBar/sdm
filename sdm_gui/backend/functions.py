# -*- coding:utf-8 -*-
import os
import tempfile
from sdm import config


TMP_DIR = config.PATH['tmp']


def download_s3(rfp, profile=None):
    suffix = rfp.split('/')[-1]
    fd, lfp = tempfile.mkstemp(suffix='.' + suffix, prefix='tmp', dir=TMP_DIR)
    cmd = ['aws', 's3', 'cp']
    if profile is not None:
        cmd += ['--profile', profile]
    cmd += [rfp, lfp, '>>{}/log'.format(TMP_DIR)]
    os.system(' '.join(cmd))
    return lfp
