# -*- coding:utf-8 -*-
from flask import Flask
from .blueprint import api_bp

app = Flask(__name__)
app.register_blueprint(api_bp, url_prefix='/api')
