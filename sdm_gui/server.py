# -*- coding: utf-8 -*-
import wx
import json
import socket
import socketserver
from threading import Thread
from . import config
from . import logger as root_logger
logger = root_logger.getChild('socketserver')


PORT_START = config.NETWORK['forntend_port_start']
PORT_END = config.NETWORK['forntend_port_end']


def core(app, message, addr):
    data = json.loads(message)
    logger.info(f"Recv: {len(message)} from {addr}")
    if app is not None:
        wx.CallAfter(
            app.canvas.structLoad,
            filename=data['type'],
            string=data['string'],
            source=f"{data['origin']} @ {addr}"
        )


class MyServerHandler(socketserver.BaseRequestHandler):
    def recvall(self, bufsize=1024):
        buf = self.request.recv(bufsize).decode()
        if buf[:6] == 'LENGTH':
            remain = int(buf.split(' ', 1)[1])
            self.request.send(buf.encode())
            raw = bytes()
            while remain > 0:
                raw += self.request.recv(bufsize)
                remain -= bufsize
            data = raw.decode()
        else:
            data = buf
        return data

    def sendall(self, msg):
        data = msg.rstrip('\n') + '\n\n'
        self.request.sendall(data.encode())

    def handle(self):
        # self.request  self.client_address  self.server
        self.sendall('SDM 200 OK')
        data = self.recvall()
        try:
            core(self.server.app, data, self.client_address)
            self.sendall('SDM 202 Accepted')
        except KeyboardInterrupt:
            raise
        except Exception as exc:
            logger.info(repr(exc))
            self.sendall('SDM 406 Not Acceptable\n' + repr(exc))
        finally:
            self.request.close()


class MyServerTask:
    def __init__(self, app):
        self.app = app
        self.port = None
        self.thread = None
        self.server = None

    def run_server(self, addr):
        port = PORT_START
        while port < PORT_END:
            try:
                self.server = socketserver.ThreadingTCPServer((addr, port), MyServerHandler)
                break
            except OSError:
                pass
            port += 1
        else:
            logger.error('Cannot start socket server')
            return
        self.port = port
        logger.info(f"listening {addr}:{port}")
        config.NETWORK['forntend_port'] = port
        config.save()
        self.server.app = self.app
        self.server.serve_forever()

    def stop_server(self):
        self.server.shutdown()
        self.server.server_close()
        config.NETWORK['forntend_port'] = 0
        config.save()

    def run_server_thread(self):
        self.thread = Thread(target=self.run_server, args=('localhost', ))
        self.thread.start()

    def stop_server_thread(self):
        self.stop_server()
        self.thread.join()


class MyClient:
    def __init__(self, host=None, port=None):
        """用来从 python 发送结构到可视化窗口服务器的接口"""
        config.load()
        self.host = host if host is not None else "127.0.0.1"
        self.port = port if port is not None else config.NETWORK['forntend_port']
        if self.port == 0:
            raise OSError("no running SDM-GUI server")
        self.client = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.client.connect((self.host, self.port))

    def recv(self):
        response = self.client.recv(1024).decode()
        rsp = response.split('\n', 1)
        head = rsp[0]
        body = rsp[1] if len(rsp) > 1 else '--'
        protocol, retcode, status = head.split(None, 2)
        assert protocol == 'SDM', 'Invalid response protocol'
        return retcode, status, body

    def sendall(self, data):
        try:
            retcode, status, body = self.recv()
            assert retcode == '200', 'Server not ready'

            raw = data.encode()
            head = f'LENGTH {len(raw)}'.encode()
            self.client.send(head)
            self.client.recv(1024)
            self.client.sendall(raw)
            retcode, status, body = self.recv()
            body = body.strip()
            logger.info(f"retcode={retcode}, status={status}, body={body}")
        finally:
            self.client.close()
