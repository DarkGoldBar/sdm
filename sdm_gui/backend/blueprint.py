# -*- coding:utf-8 -*-
from flask import request, Blueprint
from .functions import download_s3

api_bp = Blueprint('api', __name__)


@api_bp.route('/echo', methods=['POST', 'GET'])
def hello():
    if request.method == 'GET':
        response = 'hello'
    else:
        response = {
            'args': request.args,
            'form': request.form,
            'values': request.values,
            'json': request.json,
            'data': str(request.data),
        }
    return response


@api_bp.route('/v1', methods=['POST'])
def v1():
    action = request.args.get('action')
    if action == 'download':
        rfp = request.form.get('rfp', '')
        profile = request.form.get('profile', None)
        if rfp.startswith('s3://'):
            lfp = download_s3(rfp, profile)
            with open(lfp) as f:
                content = f.read()
            response = {'content': content}
    else:
        response = {'error': 'action is None'}
    return response
