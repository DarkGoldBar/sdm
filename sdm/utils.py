# -*- coding:utf-8 -*-
# Author : bochen.li
# Date : 2020/9/24
import numpy as np
import warnings


class Units:
    """
    距离基准 (Angstrom)
    能量基准 (kJ_mol)
    """
    # 常数
    Avogadro = 6.0221409e+23
    # 距离
    angstrom = 1
    bohr = 1.8897259886
    nm = 10
    m = 1e10
    # 能量
    kJ_mol = 1
    kcal_mol = 4.184
    eV = 96.4869
    hartree = 2625.5
    J = Avogadro * 1e3


def latt62Matrix(latt6):
    a = latt6[0]
    b = latt6[1]
    c = latt6[2]
    A = np.deg2rad(latt6[3])
    B = np.deg2rad(latt6[4])
    C = np.deg2rad(latt6[5])
    cos_A = np.cos(A)
    cos_B = np.cos(B)
    cos_C = np.cos(C)
    sin_C = np.sin(C)
    v_ = 1 - cos_A**2 - cos_B**2 - cos_C**2 + \
        2 * cos_A * cos_B * cos_C
    assert v_ > 0, "lattice error: alpha + beta < gamma"
    return np.array([[a, 0., 0.],
                     [b * cos_C, b * sin_C, 0.],
                     [c * cos_B,
                      c * (cos_A - cos_B * cos_C) / sin_C,
                      c * np.sqrt(v_) / sin_C]], dtype=np.float64)


def matrix2Latt6(matrix):
    a = np.linalg.norm(matrix[0])
    b = np.linalg.norm(matrix[1])
    c = np.linalg.norm(matrix[2])
    n_mat = matrix / np.array([[a, b, c]]).T
    A = np.rad2deg(np.arccos(np.dot(n_mat[1], n_mat[2])))
    B = np.rad2deg(np.arccos(np.dot(n_mat[0], n_mat[2])))
    C = np.rad2deg(np.arccos(np.dot(n_mat[0], n_mat[1])))
    return a, b, c, A, B, C


def volume_difference(struct, start=-3, stop=4, step=1):
    """体积差分
    start, stop, step: 差分起点,终点,步长 (percentage)
    """
    output = []
    for i in np.arange(start, stop, step):
        mat = struct.lattice.matrix * ((1 + i / 100.) ** (1 / 3.))
        st = struct.copy()
        st.setLattice(mat, fix_frac=True)
        output.append(st)
    return output


def angle_difference(struct, start=-3, stop=4, step=1, axis=1):
    """角度差分
    start, stop, step: 差分起点,终点,步长 (degree)
    axis: 差分角 0:alpha 1:beta 2:gamma
    """
    output = []
    for i in np.arange(start, stop, step):
        latt6 = list(struct.lattice.parameters)
        latt6[axis + 3] += i
        mat = latt62Matrix(latt6)
        st = struct.copy()
        st.setLattice(mat, fix_frac=True)
        output.append(st)
    return output


def rotation_matrix(axis, angle):
    '''Euler-Rodrigues formula for rotation matrix'''
    # Normalize the axis
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(angle / 2)
    b, c, d = -axis * np.sin(angle / 2)
    return np.array([[a*a + b*b - c*c - d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                        [2*(b*c+a*d), a*a + c*c - b*b - d*d,  2*(c*d-a*b)],
                        [2*(b*d-a*c), 2*(c*d+a*b), a*a + d*d - b*b - c*c]])


def rotate(coords, aix_a, aix_b, theta):
    """将所有点绕轴转指定角度

    Args:
        coords (numpy.array): 待旋转的点
        aix_a (numpy.array): 轴起始点
        aix_b (numpy.array): 轴终止点
        theta (float): 旋转角 (radians)

    Returns:
        numpy.array: 转后的点
    """
    _coords = np.ones((coords.shape[0], 4))
    _coords[:, :-1] = coords
    if aix_a is None:
        aix_a = np.zeros(3)

    aix = aix_b - aix_a
    aix = aix / np.linalg.norm(aix)
    n1, n2, n3 = aix.tolist()
    n12, n22, n32 = n1**2, n2**2, n3**2
    cosq, sinq = np.cos(theta), np.sin(theta)

    vec4 = np.zeros(4)
    vec4[:3] = aix_a

    t1 = np.array([
        [      n12 + (1 - n12)*cosq, n1*n2*(1 - cosq) + n3*sinq, n1*n3*(1 - cosq) - n2*sinq, 0],
        [n1*n2*(1 - cosq) - n3*sinq,       n22 + (1 - n22)*cosq, n2*n3*(1 - cosq) + n1*sinq, 0],
        [n1*n3*(1 - cosq) + n2*sinq, n2*n3*(1 - cosq) - n1*sinq,       n32 + (1 - n32)*cosq, 0],
        [                         0,                          0,                          0, 1],
    ])

    _coords -= vec4
    _coords = np.dot(_coords, t1)
    _coords += vec4

    if len(coords.shape) == 1:
        return _coords[0, :-1]
    else:
        return _coords[:, :-1]


def frac_distance2(fc1, fc2, matrix):
    """计算周期性边界下的最短距离, 仅限小距离快速计算

    Args:
        fc1, fc2: 两组分数坐标 (n*3, m*3)
        matrix: 周期性边界矩阵 (3*3)

    Return:
        dist2: 距离平方的矩阵 (n*m)
    """
    ffc1 = fc1 - np.floor(fc1)
    ffc2 = fc2 - np.floor(fc2)
    frac_vec = np.abs(ffc1.reshape(-1, 1, 3) - ffc2.reshape(1, -1, 3))
    frac_vec = 0.5 - np.abs(frac_vec - 0.5)
    cart_vec = frac_vec @ matrix
    dist2 = np.sum(cart_vec**2, axis=2)
    return dist2


def atom_idx_to_title(atom, idx):
    if idx >= 1296:
        warnings.warn("atom_idx_to_title: atoms overload")
    symbol = (atom.symbol + "_")[:2]
    base36 = "{:>2s}".format(np.base_repr(idx, 36)).replace(' ', '0').lower()
    return symbol + base36
