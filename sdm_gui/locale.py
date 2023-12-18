# -*- coding: utf-8 -*-
import os
import yaml
from . import config
from . import logger as root_logger
logger = root_logger.getChild('canvas')

# Shortcut Example
# self.myMenuItem.SetItemLabel("My &item\tCTRL+I")


def translate_windows(main, lang='cn'):
    fp = config.PATH['module'] / 'resource' / 'locale_{}.yaml'.format(lang)
    if not os.path.exists(fp):
        logger.error("locale file not found!")
        return
    with open(fp, encoding="utf-8") as f:
        loc_yaml = yaml.safe_load(f)
    for line in loc_yaml["translate"]:
        if 100 < line["id"] < 300:
            item = main.menubar.FindItemById(line["id"])
            if item is None:  # Menu 根目录对象没有ID
                pass
            else:
                item.SetItemLabel(line["label"])
                if "help" in line:
                    item.SetHelp(line["help"])
