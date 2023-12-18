# -*- coding:utf-8 -*-
from . import app
from .. import config


app.run(
    host=config.NETWORK['backend_host'],
    port=config.NETWORK['backend_port'],
)
