import os
import wx
from .locale import translate_windows
from .layout import MyApp
from .main import IMyMainFrame
from .server import MyServerTask
from . import config
from . import logger as root_logger
logger = root_logger.getChild('App')


class IMyApp(MyApp):
    def __init__(self, *args):
        self.server = None
        MyApp.__init__(self, *args)

    def InitLocale(self):
        if os.name == 'nt':
            # Walk around WinNT locale bug with pandas
            self.ResetLocale()

    def OnInit(self):
        dependence_check()
        logger.info("PATH of sdm: {}".format(config.PATH['module']))
        pos = config.PREFER['position']
        size = config.PREFER['size']
        self.main = IMyMainFrame(None, wx.ID_ANY, "sdm", pos=pos, size=size)
        self.main.SetSize(*size)
        self.main.Bind(wx.EVT_CLOSE, self.OnExitApp)
        self.SetTopWindow(self.main)
        translate_windows(main=self.main, lang="cn")
        self.main.Show()
        self.canvas = self.main.canvas
        self.server = MyServerTask(self)
        self.server.run_server_thread()
        return True

    def OnExitApp(self, evt):
        if self.server is not None:
            self.server.stop_server_thread()
            logger.info('sdm socketserver shutdown')
        pos = self.main.GetPosition()
        size = self.main.GetSize()
        config.PREFER['position'] = [pos.x, pos.y]
        config.PREFER['size'] = [size.x, size.y]
        config.save()
        self.main.Destroy()


def dependence_check():
    import sys, OpenGL  # noqa=E104
    logger.info("Dependence version check ...")
    logger.info("Python: {}".format(sys.version.replace("\n", ";")))
    logger.info("PyOpenGL: {}".format(OpenGL.__version__))
    logger.info("wxPython: {}".format(wx.__version__))


def main(init_struct=None):
    sdm = IMyApp(0)
    sdm.MainLoop()


if __name__ == "__main__":
    main()
