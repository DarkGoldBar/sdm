# -*- coding: utf-8 -*-
import requests
from .. import config


BACKEND_ADDR = config.NETWORK['backend_host'] + ':' + str(config.NETWORK['backend_port'])


def get_s3_file_string(rfp, profile, addr=BACKEND_ADDR):
    url = 'http://' + addr + '/api/v1?action=download'
    resp = requests.post(url, data={'rfp': rfp, 'profile': profile})
    if resp.status_code == 200:
        data = resp.json()
        return data['content']
