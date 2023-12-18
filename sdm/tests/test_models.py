# -*- coding:utf-8 -*-
from pathlib import Path
import numpy as np  # noqa
import sdm
from sdm.models import SpaceGroup, SymmOp  # noqa

fp = Path(__file__).parent / 'ACSALA.cif'


def test(cmds):
    mol = sdm.open(fp, astype='mol')  # noqa
    cry = sdm.open(fp)  # noqa
    err = []
    for i, cmd in enumerate(cmds):
        print('{}: {}'.format(i, cmd))
        try:
            print(eval(cmd))
        except KeyboardInterrupt:
            raise
        except Exception as e:
            print('\033[91m\033[1m' + repr(e) + '\033[0m')
            err.append((cmd, repr(e)))
        print()
    print('========== Summery ==========')
    for k, v in err:
        print(k, v)


TESTCMD1 = [
    "mol.atoms",
    "mol.bonds",
    "mol.atom_index",
    "mol.atom_titles",
    "mol.elements",
    "mol.natom",
    "mol.mass",
    "mol.cycles",
    "mol.groups",
    "mol.getConnectedGroups()",
    "mol.getNeighbor([0])",
    "mol.getBonds([[1, 2]])",
    "mol.findAtomByTitle('C1')",
    "mol.findRotatableBonds()",
    "mol.isRotatableBond(0, 1)",
    "mol.isChiralAtom(0)",
    "mol.addAtoms(['Cl'], coords=[[0, 0, 0]], charge=[1])",
    "mol.addAtomGraph(mol)",
    "mol.delAtoms([22, 23, 24, 25])",
    "mol.getAtoms([1,2,3,4,5])",
    "mol.copy()",
]

TESTCMD2 = [
    "mol.getDistance(1, 2)",
    "mol.getAngle(1, 2, 3)",
    "mol.getDihedral(1, 2, 3, 4)",
    "mol.getCentroid()",
    "mol.getGeometryCenter()",
    "mol.translate([0,1,0])",
    "mol.rotateWithAxisAngle(np.zeros(3), np.ones(3), 0.5)",
    "mol.rotateWithMatrix(np.eye(3)[[0,2,1]])",
]

TESTCMD3 = [
    "cry.frac_coords",
    "cry.volume",
    "cry.moveMolCenterInCell()",
    "cry.moveAtomsInCell()",
    "cry.moveAtomsByMolecule()",
    "cry.setLatticeLowerTriangular()",
]

TESTCMD4 = [
    "cry.lattice.copy()",
    "cry.lattice.fromParameters(4, 5, 6, 30, 40, 50)",
    "cry.lattice.parameters",
    "cry.lattice.abc",
    "cry.lattice.length",
    "cry.lattice.angles",
    "cry.lattice.inv_matrix",
    "cry.lattice.reciprocal_lattice",
    "cry.lattice.reciprocal_lattice_crystallographic",
    "cry.lattice.volume",
    "cry.lattice.height",
    "cry.lattice.metric_tensor",
    "cry.lattice.getCartesianCoords(np.eye(3) * [1, 2, 3])",
    "cry.lattice.getFractionalCoords(np.eye(3) * [0.1, 0.2, 0.3])",
    "cry.lattice.getSupercellPacking(15, 8000)",
    "cry.lattice.getLatticeSystem(14)",
    "cry.lattice.getLatticeSystemByParameters()",
    "cry.lattice.setParametersAuto(14)",
    "cry.lattice.setParameters([4, 5, 6, 30, 40, 50])",
    "cry.lattice.d_hkl([1, 2, 1])",
]


TESTCMD5 = [
    "cry.getAtoms([1, 2, 3])",
    "cry.copy()",
    "cry.getP1Structure()",
    "cry.getAsymmetricStructure()",
    "cry.getPrimitiveStructure()",
    "cry.getSupercell([3, 3, 3])",
    "cry.findEquivalentAtomAndSymmetry()",
    "cry._getVaspStyleMap()",
    "cry.delDuplicatedAtoms()",
    "cry.transformSpaceGroup(SpaceGroup.fromHallNumber(82))",
]


TESTCMD6 = [
    "cry.space_group.int_symbol",
    "cry.space_group.data",
    "cry.space_group.int_number",
    "cry.space_group.hall_number",
    "cry.space_group.int_symbol_full",
    "cry.space_group.universal_symbol",
    "cry.space_group.hall_symbol",
    "cry.space_group.lattice_system",
    "cry.space_group.copy()",
    "cry.space_group.transformCentering('C')",
    "cry.space_group.getCartesianOps(np.eye(3)*10)",
    "SpaceGroup.fromCartesianOps(np.eye(3)*10, [np.eye(4)])",
    "SpaceGroup.fromXyzOps(['x,y,z', '-x,-y,-z'])",
    "SpaceGroup.fromIntSymbol('P1c1')",
    "SpaceGroup.fromHallNumber(529)",
    "SpaceGroup.fromIntNumber(227)",
    "SpaceGroup.fromSymbol('P121')",
    "SpaceGroup.fromSymbol('P121').getTransformationMatrixAndVector('P112')",
    "cry.space_group.symmetry_ops[0].rotation_matrix",
    "cry.space_group.symmetry_ops[0].translation_vector",
    "cry.space_group.symmetry_ops[0].inverse",
    "cry.space_group.symmetry_ops[0].sign",
    "cry.space_group.symmetry_ops[0].operate([1, 2, 3])",
    "cry.space_group.symmetry_ops[0].operate_multi(np.eye(3))",
    "cry.space_group.symmetry_ops[0].as_xyz_string",
    "SymmOp.from_xyz_string('-x,-y+1/2,x+z')",
]


TESTCMD7 = [
    "cry.toString('cif')",
    "cry.toString('res')",
    "cry.toString('mol')",
    "cry.toString('xyz')",
    "cry.toString('gen')",
    "cry.toString('poscar')",
    "cry.toJson()"
]


test(TESTCMD1 + TESTCMD2 + TESTCMD3 + TESTCMD4 + TESTCMD5 + TESTCMD6 + TESTCMD7)
