# -*- coding:utf-8 -*-
import os
import yaml
import tempfile
from pathlib import Path
from . import logger as root_logger
logger = root_logger.getChild('Config')


CONFIG_FILE = 'sdm.conf'

_default_PREFER = {
    'last': '',
    'size': (900, 900),
    'position': (100, 100),
}

_default_NETWORK = {
    'backend_host': '0.0.0.0',
    'backend_port': 7701,
    'forntend_host': '127.0.0.1',
    'forntend_port': 0,
    'forntend_port_start': 7707,
    'forntend_port_end': 7727,
}


def get_path() -> dict:
    conf_dir = temp_dir = tempfile.gettempdir()
    if 'USERPROFILE' in os.environ:  # winNT
        conf_dir = Path(os.environ['USERPROFILE'])
    elif 'HOME' in os.environ:  # Linux
        conf_dir = Path(os.environ['HOME']) / '.config'
    if 'SDM_TEMP' in os.environ:
        temp_dir = os.environ['SDM_TEMP']
    if not os.path.isdir(conf_dir):
        os.mkdir(conf_dir)
    conf_file = conf_dir / CONFIG_FILE
    module_dir = Path(__file__).parent
    return {
        'config': conf_file,
        'module': module_dir,
        'tmp': temp_dir
    }


class Config():
    _save_keys = ["PREFER", "NETWORK"]

    def __init__(self):
        self.PATH = get_path()
        self.PREFER = _default_PREFER.copy()
        self.NETWORK = _default_NETWORK.copy()

    def __getitem__(self, key):
        return self.__getattribute__(key)

    def load(self):
        config = {}
        if os.path.exists(self.PATH['config']):
            with open(self.PATH['config'], encoding='utf-8') as f:
                config_str = f.read()
            config = yaml.safe_load(config_str)
        else:
            self.save()
        for k in self._save_keys:
            self[k].update(config.get(k, {}))

    def save(self):
        config = {k: self[k] for k in self._save_keys}
        with open(self.PATH['config'], 'w') as f:
            f.write(yaml.safe_dump(config))

    @classmethod
    def get_config(cls):
        cfg = cls()
        logger.info('config file: {}'.format(cfg.PATH['config']))
        try:
            cfg.load()
        except Exception:
            logger.error('config file broken, try remove sdm.conf')
        return cfg

config = Config.get_config()  # noqa
