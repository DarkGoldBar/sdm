# -*- coding: utf-8 -*-
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.arrays import vbo
from OpenGL.GL import shaders
import re
import numpy as np
from sdm.utils import rotation_matrix


class MyVBO():
    """OpenGL 顶点缓冲控制"""
    def __init__(self):
        self.vbo = None
        self.array = None
        self.length = 0
        self._enabled = False

    def update(self, data=None):
        if data is None:
            # disable current vbo
            self._enabled = False
            self.array = None
        elif len(data) == 0:
            # disable current vbo
            self._enabled = False
            self.array = None
        else:
            # set/update self.vbo
            if self.vbo is None:
                self.vbo = vbo.VBO(data)
            else:
                self.vbo.set_array(data)
            self.array = data
            self.length = data.shape[0]
            self._enabled = True

    def bind(self):
        return self.vbo.bind()

    def unbind(self):
        return self.vbo.unbind()

    def __del__(self):
        del self.vbo
        self.vbo = None


class GLShader:
    """OpenGL 着色器控制"""
    def __init__(self, vertex_code, fragment_code) -> None:
        self.program = None
        self.vertex_code = vertex_code
        self.fragment_code = fragment_code
        self.uniform = {}
        self.attrib = {}

    def Compile(self):
        # Compile
        self.program = shaders.compileProgram(
            shaders.compileShader(self.vertex_code, GL_VERTEX_SHADER),
            shaders.compileShader(self.fragment_code, GL_FRAGMENT_SHADER)
        )
        # parse vertex shader define
        ma_define = re.search(r'(.*)void main', self.vertex_code, re.DOTALL)
        attribute_key = re.findall(r'(?:in |attribute ).*\s+(\w+);', ma_define.group(1))
        uniform_key = re.findall(r'(?:uniform ).*\s+(\w+);', ma_define.group(1))
        # parse fragment shader define
        ma_define = re.search(r'(.*)void main', self.fragment_code, re.DOTALL)
        uniform_key += re.findall(r'(?:uniform ).*\s+(\w+);', ma_define.group(1))
        # allocate memory
        self.attrib = {key: glGetAttribLocation(self.program, key) for key in attribute_key}
        self.uniform = {key: glGetUniformLocation(self.program, key) for key in uniform_key}

    def __del__(self):
        if self.program is not None:
            del self.program
            self.program = None


class GLCore:
    """OpenGL 渲染核心"""
    def __init__(self):
        self.shader = None
        self.shader_selection = None
        self.matrix = {
            'ModelView': np.eye(4, dtype='f'),
            'Normal': np.eye(4, dtype='f'),
            'Projection': np.eye(4, dtype='f'),
        }
        self.vbos = {
            "triangles": MyVBO(),
            "circles": MyVBO(),
            "lines": MyVBO(),
        }

        self.is_perspective = False
        self.scale = 1.0
        self.win_wh = (640, 480)  # 宽, 高
        self.view = np.array([-10, 10, -10, 10, -500, 500], dtype=float)  # 视场左、右、上、下、前、后
        self.camera = np.array([
            [0.0, 0.0, 10.0],  # 摄像机位置
            [0.0, 0.0, 0.0],   # 摄像机目标
            [0.0, 1.0, 0.0]    # 摄像机上方
        ])

    def InitGL(self):
        glClearColor(0.0, 0.0, 0.0, 1.0)  # 设置画布背景色 RGBA
        glEnable(GL_DEPTH_TEST)
        # glEnable(GL_LINE_SMOOTH)  # 无GPU加速时很卡
        glDepthFunc(GL_LEQUAL)
        self.shader = GLShader(GLSL_vertex_shader, GLSL_fragment_shader)
        self.shader.Compile()
        self.shader_selection = GLShader(GLSL_vertex_shader_selection, GLSL_fragment_shader_selection)
        self.shader_selection.Compile()

    def Draw(self):
        # arrays -> VBO
        if self.vs is None:
            arrays = dict()
        else:
            arrays = self.vs.arrays
        for k, v in self.vbos.items():
            v.update(arrays.get(k))

        # 绘图初始化
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # 设置透视
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        wh = self.win_wh[0] / self.win_wh[1]
        if wh > 1:
            view = np.array(self.view) * [wh, wh, 1, 1, 1, 1]
        else:
            view = np.array(self.view) * [1, 1, 1/wh, 1/wh, 1, 1]
        view[:4] /= self.scale
        if self.is_perspective:
            glFrustum(*view)
        else:
            glOrtho(*view)

        # 设置摄像机
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        gluLookAt(*self.camera.flatten())
        glViewport(0, 0, *self.win_wh)

        # 读取矩阵
        glGetFloat(GL_MODELVIEW_MATRIX, self.matrix['ModelView'])
        glGetFloat(GL_PROJECTION_MATRIX, self.matrix['Projection'])
        self.matrix['Normal'] = np.linalg.inv(self.matrix['ModelView']).T

        # 着色器
        try:
            # 渲染三角形
            shader = self.shader
            buffer = self.vbos['triangles']
            if buffer._enabled:
                glUseProgram(shader.program)
                buffer.bind()
                try:
                    glUniform1i(shader.uniform['lightMode'], 1)
                    glUniform1f(shader.uniform['Ka'], 0.2)
                    glUniform1f(shader.uniform['Kd'], 0.8)
                    glUniform1f(shader.uniform['Ks'], 0.2)
                    glUniform1f(shader.uniform['shininessVal'], 1.0)
                    glUniform3f(shader.uniform['lightPos'], 200, 200, 500)
                    glUniform3f(shader.uniform['ambientColor'], 1., 1., 1.)
                    glUniform3f(shader.uniform['diffuseColor'], 1., 1., 1.)
                    glUniform3f(shader.uniform['specularColor'], 1., 1., 1.)
                    glUniformMatrix4fv(shader.uniform['mModelView'], 1, False, self.matrix['ModelView'])
                    glUniformMatrix4fv(shader.uniform['mProjection'], 1, False, self.matrix['Projection'])
                    glUniformMatrix4fv(shader.uniform['mNormal'], 1, False, self.matrix['Normal'])
                    stride = 9 * 4
                    glEnableVertexAttribArray(shader.attrib['position'])
                    glEnableVertexAttribArray(shader.attrib['normal'])
                    glEnableVertexAttribArray(shader.attrib['color'])
                    glVertexAttribPointer(shader.attrib['position'], 3, GL_FLOAT, False, stride, buffer.vbo)
                    glVertexAttribPointer(shader.attrib['normal'], 3, GL_FLOAT, False, stride, buffer.vbo + 12)
                    glVertexAttribPointer(shader.attrib['color'], 3, GL_FLOAT, False, stride, buffer.vbo + 24)
                    glDrawArrays(GL_TRIANGLES, 0, buffer.length)
                finally:
                    buffer.unbind()
                    glDisableVertexAttribArray(shader.attrib['position'])
                    glDisableVertexAttribArray(shader.attrib['normal'])
                    glDisableVertexAttribArray(shader.attrib['color'])

            # 渲染线段
            shader = self.shader
            buffer = self.vbos['lines']
            if buffer._enabled:
                glUseProgram(shader.program)
                buffer.bind()
                try:
                    glUniform1i(shader.uniform['lightMode'], 2)
                    glUniform1f(shader.uniform['Ka'], 0.8)
                    glUniform1f(shader.uniform['Kd'], 0.0)
                    glUniform1f(shader.uniform['Ks'], 0.0)
                    glUniform1f(shader.uniform['shininessVal'], 1.0)
                    glUniform3f(shader.uniform['lightPos'], 0., 0., 0.)
                    glUniform3f(shader.uniform['ambientColor'], 1., 1., 1.)
                    glUniform3f(shader.uniform['diffuseColor'], 1., 1., 1.)
                    glUniform3f(shader.uniform['specularColor'], 1., 1., 1.)
                    glUniformMatrix4fv(shader.uniform['mModelView'], 1, False, self.matrix['ModelView'])
                    glUniformMatrix4fv(shader.uniform['mProjection'], 1, False, self.matrix['Projection'])
                    glUniformMatrix4fv(shader.uniform['mNormal'], 1, False, self.matrix['Normal'])
                    stride = 6 * 4
                    glEnableVertexAttribArray(shader.attrib['position'])
                    glEnableVertexAttribArray(shader.attrib['color'])
                    glVertexAttribPointer(shader.attrib['position'], 3, GL_FLOAT, False, stride, buffer.vbo)
                    glVertexAttribPointer(shader.attrib['color'], 3, GL_FLOAT, False, stride, buffer.vbo + 12)
                    glDrawArrays(GL_LINES, 0, buffer.length)
                finally:
                    buffer.unbind()
                    glDisableVertexAttribArray(shader.attrib['position'])
                    glDisableVertexAttribArray(shader.attrib['color'])

            # 渲染选择区域
            shader = self.shader_selection
            buffer = self.vbos['circles']
            if buffer._enabled:
                glUseProgram(shader.program)
                buffer.bind()
                try:
                    glUniform4f(shader.uniform['selectColor'], 1.0, 1.0, 0.0, 0.5)
                    glUniformMatrix4fv(shader.uniform['mModelView'], 1, False, self.matrix['ModelView'])
                    glUniformMatrix4fv(shader.uniform['mProjection'], 1, False, self.matrix['Projection'])
                    stride = 6 * 4
                    glEnableVertexAttribArray(shader.attrib['position'])
                    glEnableVertexAttribArray(shader.attrib['vertexOffset'])
                    glEnableVertexAttribArray(shader.attrib['selectRadius'])
                    glVertexAttribPointer(shader.attrib['position'], 3, GL_FLOAT, False, stride, buffer.vbo)
                    glVertexAttribPointer(shader.attrib['vertexOffset'], 2, GL_FLOAT, False, stride, buffer.vbo + 12)
                    glVertexAttribPointer(shader.attrib['selectRadius'], 1, GL_FLOAT, False, stride, buffer.vbo + 20)
                    glDrawArrays(GL_TRIANGLES, 0, buffer.length)
                finally:
                    buffer.unbind()
                    glDisableVertexAttribArray(shader.attrib['position'])
                    glDisableVertexAttribArray(shader.attrib['vertexOffset'])
                    glDisableVertexAttribArray(shader.attrib['selectRadius'])
        finally:
            glUseProgram(0)
        # glutSwapBuffers()

    def glClearColor(self, *rgba):
        glClearColor(*rgba)

    def ResetCamera(self):
        self.camera = [
            [0, 0, 1],
            [0, 0, 0],
            [0, 1, 0]
        ]
        self.view[:4] = np.array([-1, 1, -1, 1], float) * self.vs.border_size
        self.scale = 1.0

    def RotateCamera(self, dx: float, dy: float):
        vec_camera = self.camera[1] - self.camera[0]
        ax_horizon = np.cross(self.camera[2], vec_camera)
        ax_vertical = np.cross(vec_camera, ax_horizon)
        ang_horizon = dy / self.win_wh[0] * 3
        ang_vertical = dx / self.win_wh[1] * 3
        mat_horizon = rotation_matrix(ax_horizon, ang_horizon)
        mat_vertical = rotation_matrix(ax_vertical, ang_vertical)
        vec_camera2 = np.dot(np.dot(vec_camera, mat_horizon), mat_vertical)

        self.camera[0] = self.camera[1] - vec_camera2
        self.camera[2] = ax_vertical / np.linalg.norm(ax_vertical)
        # glutPostRedisplay()

    def TranslateCamera(self, dx: float, dy: float):
        dist_horizon = dx / self.win_wh[0] * (self.view[1] - self.view[0])
        dist_vertical = dy / self.win_wh[1] * (self.view[1] - self.view[0])
        self.view[:4] += [dist_horizon, dist_horizon, dist_vertical, dist_vertical]
        # glutPostRedisplay()

    def getMouseRay(self, x: int, y: int) -> tuple:
        y = self.win_wh[1] - y
        ray1 = np.array(gluUnProject(x, y, 1.0))
        ray0 = np.array(gluUnProject(x, y, 0.0))
        return ray0, ray1


GLSL_vertex_shader = """#version 130
in vec3 position;
in vec3 normal;
in vec3 color;
out vec3 vertPos;
out vec3 normalInterp;
invariant out vec3 vertColor;
uniform mat4 mModelView;
uniform mat4 mNormal;
uniform mat4 mProjection;

void main(){
  vec4 vertPos4 = mModelView * vec4(position, 1.0);
  vertColor = color;
  vertPos = vec3(vertPos4) / vertPos4.w;
  normalInterp = vec3(mNormal * vec4(normal, 1.0));
  gl_Position = mProjection * vertPos4;
}
"""

GLSL_fragment_shader = """#version 130
precision mediump float;
in vec3 normalInterp;  // Surface normal
in vec3 vertPos;       // Vertex position
invariant in vec3 vertColor;     // Vertex color
uniform int lightMode;   // Rendering mode
uniform float Ka;   // Ambient reflection coefficient
uniform float Kd;   // Diffuse reflection coefficient
uniform float Ks;   // Specular reflection coefficient
uniform float shininessVal; // Shininess
// Material color
uniform vec3 ambientColor;
uniform vec3 diffuseColor;
uniform vec3 specularColor;
uniform vec3 lightPos; // Light position

void main() {
  vec3 N = normalize(normalInterp);
  vec3 L = normalize(lightPos - vertPos);

  // Lambert's cosine law
  float lambertian = max(dot(N, L), 0.0);
  float specular = 0.0;
  if(lambertian > 0.0) {
    vec3 R = reflect(-L, N);      // Reflected light vector
    vec3 V = normalize(-vertPos); // Vector to viewer
    // Compute the specular term
    float specAngle = max(dot(R, V), 0.0);
    specular = pow(specAngle, shininessVal);
  }
  if(lightMode == 1) gl_FragColor = vec4(vertColor * (Ka * ambientColor + Kd * lambertian * diffuseColor + Ks * specular * specularColor), 1.0);
  // only ambient
  if(lightMode == 2) gl_FragColor = vec4(Ka * ambientColor * vertColor, 1.0);
  // only diffuse
  if(lightMode == 3) gl_FragColor = vec4(Kd * lambertian * diffuseColor * vertColor, 1.0);
  // only specular
  if(lightMode == 4) gl_FragColor = vec4(Ks * specular * specularColor * vertColor, 1.0);
}
"""

GLSL_vertex_shader_selection = """#version 130
in vec3 position;
in vec2 vertexOffset;
in float selectRadius;
out vec3 selectVertPos;
invariant out float radius;
invariant out vec3 selectCenter;
uniform mat4 mModelView;
uniform mat4 mProjection;

void main(){
  vec4 pos4 = mModelView * vec4(position, 1.0);
  selectCenter = (vec3(pos4) / pos4.w);
  selectVertPos = vec3(pos4.xy + vertexOffset.xy, pos4.z);
  radius = selectRadius;
  gl_Position = mProjection * vec4(selectVertPos, pos4.w);
}
"""

GLSL_fragment_shader_selection = """#version 130
precision mediump float;
in vec3 selectVertPos;
invariant in vec3 selectCenter;
invariant in float radius;
uniform vec4 selectColor;

void main() {
  if (length(selectVertPos - selectCenter) <= (radius)) {
    gl_FragColor = selectColor;
  }
  else {
    discard;
  }
}
"""
