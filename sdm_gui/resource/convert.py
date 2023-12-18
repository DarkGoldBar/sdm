#!/home/bochen.li/.local/bin/python3
# -*- coding: utf-8 -*-

from PIL import Image
def convert32x32(fp):
    img = Image.open(fp)
    img2 = img.resize((32, 32))
    img2.save(fp.rsplit('.', 1)[0] + '_32.png')


if __name__ == '__main__':
    import sys
    for fp in sys.argv[1:]:
        convert32x32(fp)
