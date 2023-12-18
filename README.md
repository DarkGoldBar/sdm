## SDM
SDM 的全称是 San-D Molecule  

#### 主要特性
- 基于 Dict of dict 的 networkx.Graph 数据结构，更易于分子、晶体编辑
- 晶体的分数坐标(frac_coords)和直角坐标(coords)相互关联, 不需要手动更新
- 基于wxPython的可视化编辑器

#### 主要依赖
- 数据结构
    - numpy
    - networkx
- 可视化
    - wxPython
    - PyOpenGL

## 安装
#### 只使用models部分
```
python setup.py install
```
```
pip install git+ssh://git@bitbucket.org/xtalpi/sdm.git
```

#### 使用可视化
可视化已测试的系统:
conda-python3.8 @ win10
conda-python3.8\3.9 @ Ubuntu18.04
conda-python3.8 @ Ubuntu20.04
conda-python3.9 @ WSL (拖拽打开功能失效)

1. 使用pip安装
```
apt install libgtk-3-dev
pip install -r requirements.txt
```

2. 安装sdm
`python setup.py install`

3. 运行可视化界面
`python3 -m sdm_gui`



## 使用模块
#### 晶体、分子对象读写
```
import sdm
struct = sdm.open('path/to/file.cif')
mol = sdm.open('path/to/file.xyz')

struct.save('path/to/file.cif')
mol.save('path/to/file.cif')
```

#### 晶体、分子对象可视化
*在可视化界面已经打开的情况下*
`struct.show()`

## models数据结构
```
(Structure)
├── space_group (SpaceGroup)
│   └── symmops (List[SymmOp])
│
├── lattice (Lattice)
│   └── matrix (ndarray)
│
(inherit:Molecule)
│
├── coords(ndarray)
│
└── graph(nx.Graph)
    ├── nodes (dict-like)
    │   ├── Key: index (int)
    │   ├── Value: atom (Atom)
    │   │   ├── ele (Element(Enum))
    │   │   └── _private_values  # Override value in ele
    │   │
    │   ├── Value: title (str)
    │   └── Value: ... # user defined attributes
    │
    └── edges (dict-like)
        ├── Key: index, index (int)
        ├── Value: type (BondType(Enum))
        └── Value: ... # user defined attributes
```

## 可视化编辑器

#### 按键绑定
左键单击 选择原子
左键双击 选择分子
右键单击 快捷菜单
Shift + 左键 增加选择
右键旋转
中键平移

Delete 删除选中
PageUp 上一结构
PageDown 下一结构

## CHANGELOG
#### v0.2.0
增加 找键方法 idatm
增加 内坐标支持
增加 生成slab
文件支持 dmol3: .incoor 文件读写
支持打开s3文件

#### v0.3.0
分离 sdm & sdm_gui

#### v0.3.2
修复 PBE getLatticeSystem 三方晶系问题
修复 AtomGraph.getAtoms
修复 很多重构后的错误

#### TODO
框选
自动补氢
原子选中补氢
原子选中变更属性

状态栏属性变为可选择文字


超分子 SuperMolecule, 用于分子间作用分析
增加文件支持
