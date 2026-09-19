# -*- coding: utf-8 -*-
# =============================================================================
#  l12_model_builder.py
#  Stand-alone Abaqus/CAE model builder for the hinged composite-wing model
#  (original deck: l12.inp, Abaqus/CAE 2017).
#
#  This script is NOT a recorded macro and it does NOT call
#  mdb.ModelFromInputFile(). Every material, section, composite layup, set,
#  surface, constraint, interaction, step, load, BC, output request and job is
#  written out explicitly below, so you can read and edit the model in Python.
#  The only thing taken from the deck is the raw mesh table (nodal coordinates
#  and element connectivity), because the wing is a swept, twisted, non-
#  prismatic surface that exists only as a mesh - it cannot be re-created
#  parametrically. That table is read by the small reader in SECTION 1
#  (plain Python, no Abaqus import machinery).
#
#  RUN:  Abaqus/CAE -> File -> Run Script...      (GUI)
#        abaqus cae script=l12_model_builder.py   (GUI, autorun)
#        abaqus cae noGUI=l12_model_builder.py    (batch)
#
#  Known API/modeling issues found during review are corrected below.
# =============================================================================

from abaqus import *
from abaqusConstants import *
import mesh, regionToolset, section
import os, sys, math

# ============================== SECTION 0 - SETTINGS =========================
MESH_FILE       = r'l12.inp'   # source of nodes/elements only
MODEL_NAME      = 'l12_fixed'
JOB_NAME        = 'l12_fixed'
N_CPUS          = 4

# Safer GUI workflow:
#   False -> build/check the model and write the input deck, but do not submit.
#   True  -> submit automatically after the model is built.
SUBMIT_JOB      = False
WAIT_FOR_JOB    = True

FIRST_ORDER     = True   # FIX #1: drop mid-side nodes  C3D20R->C3D8R,
                         # S8R->S4R, STRI65->S3R. Required for Hashin
                         # progressive damage + element deletion + cohesive
                         # contact, which are not supported on 2nd-order
                         # shells/solids in Abaqus/Standard.
PROGRESSIVE_DAMAGE = True  # Hashin initiation + energy evolution + deletion
BOND_TYPE       = 'COHESIVE'  # 'COHESIVE' (debonding study) or 'TIE' (rigid bond)
                              # FIX #2: never both at the same time.
LOAD_SCALE      = 1.0
# loads from the original deck (N, N-mm)
TIP_FORCE       = 161.4
ROOT_MOMENT     = 254153.0

# ======================= SECTION 1 - MINIMAL MESH READER =====================
def _fl(s):
    return float(s)

class Deck(object):
    """Very small Abaqus-keyword reader: parts (nodes/elements/elsets),
    assembly instances, assembly nodes, n/elsets and element surfaces."""
    def __init__(self, path):
        self.parts = {}          # part -> {'nodes':{id:(x,y,z)}, 'elems':{type:{id:[conn]}}, 'elsets':{name:[ids]}}
        self.instances = {}      # inst -> {'part':p, 'trans':(dx,dy,dz), 'rot':(pt1,pt2,ang) or None}
        self.anodes = {}         # assembly-level reference nodes
        self.nsets = {}          # (name) -> {'inst':inst or None, 'ids':[...]}
        self.elsets = {}         # (name) -> {'inst':inst or None, 'ids':[...]}
        self.surfaces = {}       # name -> [(elsetname, face)]
        self._read(path)

    @staticmethod
    def _opts(line):
        o = {}
        for tok in line.split(',')[1:]:
            if '=' in tok:
                k, v = tok.split('=', 1)
                o[k.strip().lower()] = v.strip().strip('"')
            else:
                o[tok.strip().lower()] = True
        return o

    def _read(self, path):
        f = open(path, 'r')
        lines = f.read().split('\n')
        f.close()
        part = None; inst = None; mode = None; ctx = None
        pend = []          # continuation buffer for elements
        i = 0
        while i < len(lines):
            ln = lines[i].rstrip(); i += 1
            if not ln.strip():
                continue
            if ln.startswith('**'):
                continue
            if ln.startswith('*'):
                kw = ln.split(',')[0][1:].strip().upper()
                o = self._opts(ln)
                mode = None; pend = []
                if kw == 'PART':
                    part = o['name']; self.parts[part] = {'nodes': {}, 'elems': {}, 'elsets': {}}
                elif kw == 'END PART':
                    part = None
                elif kw == 'INSTANCE':
                    inst = o['name']
                    self.instances[inst] = {'part': o['part'], 'trans': (0., 0., 0.), 'rot': None, 'lines': []}
                elif kw == 'END INSTANCE':
                    d = self.instances[inst]['lines']
                    if len(d) >= 1:
                        v = [_fl(x) for x in d[0].split(',') if x.strip()]
                        self.instances[inst]['trans'] = tuple(v[:3])
                    if len(d) >= 2:
                        v = [_fl(x) for x in d[1].split(',') if x.strip()]
                        self.instances[inst]['rot'] = (tuple(v[0:3]), tuple(v[3:6]), v[6])
                    inst = None
                elif kw == 'NODE':
                    mode = 'node'
                elif kw == 'ELEMENT':
                    mode = 'elem'; ctx = o['type'].upper()
                    if part:
                        self.parts[part]['elems'].setdefault(ctx, {})
                elif kw == 'NSET':
                    mode = 'nset'
                    ctx = (o['nset'], o.get('instance'), 'generate' in o)
                    self.nsets.setdefault(o['nset'], {'inst': o.get('instance'), 'ids': []})
                elif kw == 'ELSET':
                    mode = 'elset'
                    ctx = (o['elset'], o.get('instance'), 'generate' in o)
                    tgt = self.parts[part]['elsets'] if part else self.elsets
                    if part:
                        tgt.setdefault(o['elset'], [])
                    else:
                        tgt.setdefault(o['elset'], {'inst': o.get('instance'), 'ids': []})
                elif kw == 'SURFACE':
                    mode = 'surf'; ctx = o['name']; self.surfaces.setdefault(ctx, [])
                continue
            # ---------------- data lines ----------------
            if inst is not None and mode is None:
                self.instances[inst]['lines'].append(ln); continue
            if mode == 'node':
                p = [x for x in ln.split(',') if x.strip()]
                nid = int(p[0]); xyz = tuple(_fl(x) for x in p[1:4])
                if part:
                    self.parts[part]['nodes'][nid] = xyz
                else:
                    self.anodes[nid] = xyz
            elif mode == 'elem' and part:
                pend.append(ln)
                if ln.rstrip().endswith(','):
                    continue
                p = [int(x) for x in (','.join(pend)).split(',') if x.strip()]
                pend = []
                self.parts[part]['elems'][ctx][p[0]] = p[1:]
            elif mode in ('nset', 'elset'):
                name, ins, gen = ctx
                vals = [x.strip() for x in ln.split(',') if x.strip()]
                if gen:
                    a, b = int(_fl(vals[0])), int(_fl(vals[1]))
                    s = int(_fl(vals[2])) if len(vals) > 2 else 1
                    ids = list(range(a, b + 1, s))
                else:
                    ids = [int(v) for v in vals if v.lstrip('-').isdigit()]
                if mode == 'nset':
                    self.nsets[name]['ids'].extend(ids)
                elif part:
                    self.parts[part]['elsets'][name].extend(ids)
                else:
                    self.elsets[name]['ids'].extend(ids)
            elif mode == 'surf':
                p = [x.strip() for x in ln.split(',')]
                if len(p) >= 2:
                    self.surfaces[ctx].append((p[0], p[1].upper()))

print('Reading mesh from %s ...' % MESH_FILE)
mesh_path = MESH_FILE if os.path.isabs(MESH_FILE) else os.path.join(os.getcwd(), MESH_FILE)
if not os.path.isfile(mesh_path):
    raise ValueError('Mesh file not found: %s' % mesh_path)
deck = Deck(mesh_path)
for p in deck.parts:
    print('   %-14s %6d nodes  %s' % (p, len(deck.parts[p]['nodes']),
          dict((t, len(e)) for t, e in deck.parts[p]['elems'].items())))

# corner-node counts for the first-order downgrade (FIX #1)
DOWN = {'C3D20R': ('C3D8R', 8), 'C3D20': ('C3D8R', 8),
        'S8R':   ('S4R', 4),   'S8R5': ('S4R', 4),
        'STRI65': ('S3R', 3),  'STRI3': ('S3R', 3)}
# Element deletion is assigned later, by wing element region.  This is
# important because WING_POINTS contains both deformable damage regions and
# rigid-body tip patches; a single global ElemType would incorrectly enable
# deletion on the rigid-body tip.
def _etype(code):
    return mesh.ElemType(elemCode=code, elemLibrary=STANDARD)

ETYPE = {
    'C3D8R': _etype(C3D8R),
    'S4R': _etype(S4R),
    'S3R': _etype(S3R),
    'C3D20R': _etype(C3D20R),
    'S8R': _etype(S8R),
    'STRI65': _etype(STRI65),
}

# ======================= SECTION 2 - MODEL AND MESH PARTS ====================
if MODEL_NAME in mdb.models.keys():
    del mdb.models[MODEL_NAME]
m = mdb.Model(name=MODEL_NAME, modelType=STANDARD_EXPLICIT)
if 'Model-1' in mdb.models.keys():
    try: del mdb.models['Model-1']
    except Exception: pass

def build_part(pname, dim):
    src = deck.parts[pname]
    nodes = []; nmap = {}
    keep = set()
    elems = []
    for etype, edict in src['elems'].items():
        tgt, nc = DOWN.get(etype, (etype, None))
        for eid, conn in edict.items():
            c = conn[:nc] if (FIRST_ORDER and nc) else conn
            keep.update(c)
            elems.append((eid, tgt if FIRST_ORDER else etype, c))
    for nid in sorted(keep):
        nmap[nid] = len(nodes)
        if nid not in src['nodes']:
            raise ValueError("Element references missing node %s in part '%s'." %
                             (nid, pname))
        nodes.append(src['nodes'][nid])
    conn_by_type = {}
    labels_by_type = {}
    for eid, t, c in elems:
        conn_by_type.setdefault(t, []).append([nmap[x] for x in c])
        labels_by_type.setdefault(t, []).append(eid)
    p = m.PartFromNodesAndElements(
        name=pname, dimensionality=dim, type=DEFORMABLE_BODY,
        nodes=nodes,
        elements=[tuple(conn_by_type[t]) for t in sorted(conn_by_type)],
        elementTypes=[ETYPE[t] for t in sorted(conn_by_type)],
        nodeLabels=tuple(sorted(keep)))
    # element sets defined inside the part (composite layup regions, sections)
    for sname, ids in src['elsets'].items():
        try:
            p.SetFromElementLabels(name=sname,
                                   elementLabels=tuple(sorted(set(ids))))
        except Exception:
            pass
    print('   part %-14s built (%d nodes, %d elements)' % (pname, len(nodes), len(elems)))
    return p

print('Building mesh parts ...')
for pn in ('HINGE_LEAF_A', 'HINGE_LEAF_B', 'Pin'):
    build_part(pn, THREE_D)
build_part('WING_POINTS', THREE_D)

# ======================== SECTION 3 - MATERIALS ==============================
# FIX #3: densities were missing on the metals (any dynamic/frequency/gravity
#         step would fail).  FIX #4: both metals were purely elastic although
#         PE/PEEQ/PEMAG were requested as output - plasticity added.
al = m.Material(name='AL-1100-H14')
al.Density(table=((2.71e-09, ), ))
al.Elastic(table=((70000.0, 0.33), ))
al.Plastic(table=((103.0, 0.0), (117.0, 0.02), (130.0, 0.10)))

st = m.Material(name='Steel-4130-QT')
st.Density(table=((7.85e-09, ), ))
st.Elastic(table=((200000.0, 0.3), ))
st.Plastic(table=((655.0, 0.0), (850.0, 0.08)))

ge = m.Material(name='Glass-Epoxy-Lamina')
ge.Density(table=((1.9e-09, ), ))
ge.Elastic(type=LAMINA, table=((38000.0, 8500.0, 0.27, 4500.0, 4000.0, 3800.0), ))
if PROGRESSIVE_DAMAGE:
    # Xt, Xc, Yt, Yc, Sl, St
    ge.HashinDamageInitiation(table=((1000.0, 700.0, 65.0, 200.0, 85.0, 65.0), ))
    # FIX #5: the original deck had a Hashin initiation criterion but NO
    # damage evolution, so nothing ever degraded and the
    # "ELEMENT DELETION = YES" section control could never delete anything.
    ge.hashinDamageInitiation.DamageEvolution(
        type=ENERGY, table=((12.0, 10.0, 1.0, 1.0), ))
    ge.hashinDamageInitiation.DamageStabilization(
        fiberTensileCoeff=5e-05, fiberCompressiveCoeff=5e-05,
        matrixTensileCoeff=5e-05, matrixCompressiveCoeff=5e-05)

# ==================== SECTION 4 - SECTIONS AND COMPOSITE LAYUPS ==============
m.HomogeneousSolidSection(name='SEC-HINGE-AL1100', material='AL-1100-H14', thickness=None)
m.HomogeneousSolidSection(name='SEC-PIN-4130QT',  material='Steel-4130-QT', thickness=None)

# NOTE:
# CompositeShellSection() does not accept a controlName argument.  In this
# mesh-built model, element deletion is controlled through ElemType
# (elemDeletion=ON/OFF) assigned to element regions.
#
# We explicitly separate deformable composite regions from the rigid-body tip
# regions.  This prevents a failed tip element from being deleted even though
# the same orphan-mesh part also contains damageable laminate elements.

PLY_T = 0.2
LAYUPS = {   # name : (stacking sequence, part element-set name)
    'CompositeLayup-Top':    ([0., 0., 0., 10., -10., -10., 10., 0., 0., 0.], 'CompositeLayup-Top-1'),
    'CompositeLayup-Bottom': ([0., 0., 0., 10., -10., -10., 10., 0., 0., 0.], 'CompositeLayup-Bottom-1'),
    'CompositeLayup-Root':   ([45., -45., 0., 90., 0., 0., 90., 0., -45., 45.], 'CompositeLayup-Root-1-4'),
    'CompositeLayup-Tip':    ([45., -45., 0., 90., 0., 0., 90., 0., -45., 45.], 'CompositeLayup-Tip-1-6'),
}
check_required_mesh_content()
wp = m.parts['WING_POINTS']
for lname, (angles, esetName) in LAYUPS.items():
    layers = []
    for k, ang in enumerate(angles):
        layers.append(section.SectionLayer(material='Glass-Epoxy-Lamina',
                                           thickness=PLY_T, orientAngle=ang,
                                           numIntPts=3, plyName='Ply-%d' % (k + 1)))
    m.CompositeShellSection(name=lname, preIntegrate=False, idealization=NO_IDEALIZATION,
                            symmetric=False, thicknessType=UNIFORM, poissonDefinition=DEFAULT,
                            thicknessModulus=None, temperature=GRADIENT, useDensity=OFF,
                            integrationRule=SIMPSON, layup=layers, layupName=lname)
    if esetName in wp.sets.keys():
        wp.SectionAssignment(region=wp.sets[esetName], sectionName=lname,
                             offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='')
# extra root/tip patches of the original deck
for extra, lname in (('CompositeLayup-Root-1-3', 'CompositeLayup-Root'),
                     ('CompositeLayup-Tip-1-5', 'CompositeLayup-Tip')):
    if extra in wp.sets.keys():
        wp.SectionAssignment(region=wp.sets[extra], sectionName=lname,
                             offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='')

# Progressive-damage element deletion:
#   - ON  for Top/Bottom/Root laminate regions
#   - OFF for Tip laminate regions that are tied to rigid-body RPs
#
# setElementType() is valid for orphan mesh elements/element sets.  We use the
# same first-order shell family already present in the part and change only the
# deletion flag.
if PROGRESSIVE_DAMAGE:
    damage_sets = []
    for _nm in ('CompositeLayup-Top-1',
                'CompositeLayup-Bottom-1',
                'CompositeLayup-Root-1-4',
                'CompositeLayup-Root-1-3'):
        if _nm in wp.sets.keys():
            damage_sets.append(wp.sets[_nm])

    tip_sets = []
    for _nm in ('CompositeLayup-Tip-1-6',
                'CompositeLayup-Tip-1-5'):
        if _nm in wp.sets.keys():
            tip_sets.append(wp.sets[_nm])

    # Determine which shell families are actually present in each region.
    # S4R and S3R are both supported by the downgraded mesh.
    for _reg in damage_sets:
        wp.setElementType(
            regions=(_reg,),
            elemTypes=(mesh.ElemType(elemCode=S4R,
                                     elemLibrary=STANDARD,
                                     elemDeletion=ON),
                       mesh.ElemType(elemCode=S3R,
                                     elemLibrary=STANDARD,
                                     elemDeletion=ON)))
    for _reg in tip_sets:
        wp.setElementType(
            regions=(_reg,),
            elemTypes=(mesh.ElemType(elemCode=S4R,
                                     elemLibrary=STANDARD,
                                     elemDeletion=OFF),
                       mesh.ElemType(elemCode=S3R,
                                     elemLibrary=STANDARD,
                                     elemDeletion=OFF)))
    print('Progressive-damage element deletion assigned by wing region.')

for pn, sec in (('HINGE_LEAF_A', 'SEC-HINGE-AL1100'),
                ('HINGE_LEAF_B', 'SEC-HINGE-AL1100'),
                ('Pin', 'SEC-PIN-4130QT')):
    p = m.parts[pn]
    reg = regionToolset.Region(elements=p.elements)
    p.SectionAssignment(region=reg, sectionName=sec, offset=0.0,
                        offsetType=MIDDLE_SURFACE, offsetField='')


# =========================== MODEL VALIDATION ================================
def require_part(name):
    if name not in deck.parts:
        raise ValueError("Required part '%s' is missing from %s" % (name, MESH_FILE))

def require_instance(name):
    if name not in deck.instances:
        raise ValueError("Required instance '%s' is missing from %s" % (name, MESH_FILE))

def require_part_set(part_name, set_name):
    p = m.parts[part_name]
    if set_name not in p.sets.keys():
        raise ValueError("Required part element set '%s' is missing on '%s'." %
                         (set_name, part_name))

def require_assembly_set(name):
    if name not in a.sets.keys():
        raise ValueError("Required assembly set '%s' was not created." % name)

def check_required_mesh_content():
    for pn in ('HINGE_LEAF_A', 'HINGE_LEAF_B', 'Pin', 'WING_POINTS'):
        require_part(pn)
    for ins in ('HINGE_LEAF_A-1', 'HINGE_LEAF_B-1',
                'Pin-1', 'WING_POINTS-1', 'WING_POINTS-2'):
        require_instance(ins)

    wp_sets = set(m.parts['WING_POINTS'].sets.keys())
    for _, eset in LAYUPS.values():
        if eset not in wp_sets:
            raise ValueError("Wing layup set '%s' is missing." % eset)

    if 46521 not in deck.parts['WING_POINTS']['nodes']:
        raise ValueError("WING_POINTS reference node 46521 is missing.")

    print('Model-content validation: OK')


# ========================= SECTION 5 - ASSEMBLY ==============================
a = m.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
for iname, idat in deck.instances.items():
    inst = a.Instance(name=iname, part=m.parts[idat['part']], dependent=ON)
    t = idat['trans']
    if any(abs(v) > 1e-12 for v in t):
        a.translate(instanceList=(iname, ), vector=t)
    if idat['rot']:
        p1, p2, ang = idat['rot']
        a.rotate(instanceList=(iname, ), axisPoint=p1, axisDirection=
                 tuple(p2[k] - p1[k] for k in range(3)), angle=ang)
print('Assembly instances : %s' % ', '.join(sorted(a.instances.keys())))

# ---- reference points (they were plain assembly nodes in the deck) ----------
RP = {}
RP_DEF = {                      # name          : (x, y, z)
    'RP-PIN-CENTER'  : (0.0,   0.0, 208.5),      # node 2  - pin coupling / BCs
    'RP-HINGE-RIGHT' : (30.0, -15.0, 127.5),     # node 4  - right root
    'RP-HINGE-LEFT'  : (-30.0, -15.0, 127.5),    # node 7  - left  root
    'RP-TIP-RIGHT'   : (2454.98682, -23.112175, 102.509102),   # node 5
    'RP-TIP-LEFT'    : (-2454.97778, -4.54690933, 152.490906), # node 6
}
# FIX #7: the deck also contained two unused / duplicated reference points
#         (node 1 and "Rroot" node 3, coincident with node 4) - dropped.
def _place(iname, xyz):
    """Apply the instance translation/rotation of the deck to a part point."""
    d = deck.instances[iname]
    x = list(xyz)
    if d['rot']:
        p1, p2, ang = d['rot']
        ax = [p2[k] - p1[k] for k in range(3)]
        L = math.sqrt(sum(v * v for v in ax)) or 1.0
        u = [v / L for v in ax]
        th = math.radians(ang); c, sn = math.cos(th), math.sin(th)
        r = [x[k] - p1[k] for k in range(3)]
        dot = sum(u[k] * r[k] for k in range(3))
        cr = [u[1] * r[2] - u[2] * r[1], u[2] * r[0] - u[0] * r[2], u[0] * r[1] - u[1] * r[0]]
        x = [p1[k] + r[k] * c + cr[k] * sn + u[k] * dot * (1 - c) for k in range(3)]
    t = d['trans']
    return tuple(x[k] + t[k] for k in range(3))

# wing reference points (part node 46521 = WING_POINTS-RefPt_)
_wrp = deck.parts['WING_POINTS']['nodes'][46521]
RP_DEF['RP-WING-RIGHT'] = _place('WING_POINTS-1', _wrp)
RP_DEF['RP-WING-LEFT']  = _place('WING_POINTS-2', _wrp)

for nm, xyz in RP_DEF.items():
    f = a.ReferencePoint(point=xyz)
    RP[nm] = a.referencePoints[f.id]
    a.Set(name=nm, referencePoints=(RP[nm], ))

# ---- helper: build assembly sets / surfaces from the deck's label lists -----
def inst_elems(iname, labels):
    return a.instances[iname].elements.sequenceFromLabels(tuple(sorted(set(labels))))

def inst_nodes(iname, labels):
    return a.instances[iname].nodes.sequenceFromLabels(tuple(sorted(set(labels))))

FACEKEY = {'S1': 'face1Elements', 'S2': 'face2Elements', 'S3': 'face3Elements',
           'S4': 'face4Elements', 'S5': 'face5Elements', 'S6': 'face6Elements',
           'SPOS': 'side1Elements', 'SNEG': 'side2Elements'}

def make_surface(newname, deckSurfName):
    """Rebuild an element-based surface from the deck definition."""
    kw = {}
    for esetName, face in deck.surfaces.get(deckSurfName, []):
        ent = deck.elsets.get(esetName)
        if not ent or not ent['inst']:
            continue
        key = FACEKEY.get(face)
        if key is None:
            continue
        seq = inst_elems(ent['inst'], ent['ids'])
        kw[key] = (kw[key] + seq) if key in kw else seq
    if not kw:
        print('   ! surface %s could not be rebuilt' % deckSurfName)
        return None
    kw['name'] = newname
    return a.Surface(**kw)

SRF = {}
SRF['WING-R-BOND']  = make_surface('WING-R-BOND',  '_PickedSurf695')   # right wing root patch
SRF['HINGE-A-BOND'] = make_surface('HINGE-A-BOND', '_PickedSurf696')   # hinge leaf A face
SRF['WING-L-BOND']  = make_surface('WING-L-BOND',  '_PickedSurf697')   # left wing root patch
SRF['HINGE-B-BOND'] = make_surface('HINGE-B-BOND', '_PickedSurf698')   # hinge leaf B face
SRF['PIN-BORE']     = make_surface('PIN-BORE',     '_PickedSurf707')   # pin barrel
SRF['WING-R-RP']    = make_surface('WING-R-RP',    '_PickedSurf682')   # right wing coupling patch
SRF['WING-L-RP']    = make_surface('WING-L-RP',    '_PickedSurf680')   # left  wing coupling patch
SRF['HINGE-A-TIE']  = make_surface('HINGE-A-TIE',  '_PickedSurf729')
SRF['WING-R-ROOT']  = make_surface('WING-R-ROOT',  '_PickedSurf724')   # right root patch
SRF['WING-L-ROOT']  = make_surface('WING-L-ROOT',  '_PickedSurf726')   # left  root patch

# node sets used by BCs / rigid bodies
def aset(name, deckSetName):
    ent = deck.nsets.get(deckSetName)
    if ent and ent['inst']:
        return a.Set(name=name, nodes=inst_nodes(ent['inst'], ent['ids']))
    return None

aset('PIN-BC-NODES', '_PickedSet699')
aset('WING-R-TIP',   '_PickedSet704')
aset('WING-L-TIP',   '_PickedSet712')

for _s in ('PIN-BC-NODES', 'WING-R-TIP', 'WING-L-TIP'):
    require_assembly_set(_s)

for _s in ('RP-PIN-CENTER', 'RP-HINGE-RIGHT', 'RP-HINGE-LEFT',
           'RP-TIP-RIGHT', 'RP-TIP-LEFT', 'RP-WING-RIGHT', 'RP-WING-LEFT'):
    require_assembly_set(_s)

# ======================= SECTION 6 - CONSTRAINTS =============================
# FIX #8 (the important one).  In the original deck the 77 root-patch shell
# elements of each wing (10375-10451) were simultaneously
#     (a) slave nodes of a *Rigid Body  (R5-RightRoot / R8-"LeftWingtip"),
#     (b) the slave surface of a *Tie   (Constraint-8 / Constraint-9),
#     (c) the cohesive *Contact Pair    (Int-1 / Int-2).
# That is a triple over-constraint on the same nodes: the rigid body wins, the
# tie and the cohesive bond are silently killed (or Abaqus aborts with
# over-constraint / zero-pivot errors), so the joint could never debond - the
# whole point of the model.  Here the root patches stay deformable, load is
# introduced through DISTRIBUTING couplings, and the hinge-to-skin joint is
# made EITHER cohesive OR tied, never both.
m.Coupling(name='R5-RightRoot', controlPoint=a.sets['RP-HINGE-RIGHT'],
           surface=SRF['WING-R-ROOT'], influenceRadius=WHOLE_SURFACE,
           couplingType=DISTRIBUTING, weightingMethod=UNIFORM,
           localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
m.Coupling(name='R8-LeftRoot', controlPoint=a.sets['RP-HINGE-LEFT'],
           surface=SRF['WING-L-ROOT'], influenceRadius=WHOLE_SURFACE,
           couplingType=DISTRIBUTING, weightingMethod=UNIFORM,
           localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
# wing load-introduction couplings (kinematic, as in the original deck)
m.Coupling(name='Rp-RightWing', controlPoint=a.sets['RP-WING-RIGHT'],
           surface=SRF['WING-R-RP'], influenceRadius=WHOLE_SURFACE,
           couplingType=KINEMATIC, localCsys=None,
           u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
m.Coupling(name='Rp-LeftWing', controlPoint=a.sets['RP-WING-LEFT'],
           surface=SRF['WING-L-RP'], influenceRadius=WHOLE_SURFACE,
           couplingType=KINEMATIC, localCsys=None,
           u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
# pin -> its own reference point (kinematic, as in the original)
m.Coupling(name='R4-CenterPin', controlPoint=a.sets['RP-PIN-CENTER'],
           surface=SRF['PIN-BORE'], influenceRadius=WHOLE_SURFACE,
           couplingType=KINEMATIC, localCsys=None,
           u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
# wing tips: rigid bodies are fine there (no contact, no bond, no failure)
m.RigidBody(name='RB-TipRight', refPointRegion=a.sets['RP-TIP-RIGHT'],
            tieRegion=a.sets['WING-R-TIP'])
m.RigidBody(name='RB-TipLeft',  refPointRegion=a.sets['RP-TIP-LEFT'],
            tieRegion=a.sets['WING-L-TIP'])

# ======================= SECTION 7 - INTERACTIONS ============================
if BOND_TYPE == 'COHESIVE':
    ip = m.ContactProperty('IntProp-1')
    ip.TangentialBehavior(formulation=PENALTY, table=((1.0, ), ),
                          fraction=0.005, maximumElasticSlip=FRACTION)
    ip.NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON,
                      constraintEnforcementMethod=DEFAULT)
    ip.CohesiveBehavior(defaultPenalties=OFF, table=((10000.0, 3700.0, 3700.0), ))
    ip.Damage(initTable=((30.0, 25.0, 25.0), ), criterion=QUAD_TRACTION,
              useEvolution=ON, evolutionType=ENERGY, mixedModeType=BK,
              exponent=1.75, evolTable=((0.35, 0.9, 1.75), ),
              useStabilization=ON, viscosityCoef=5.5e-05)
    # FIX #9: slave/master were reversed (the stiff aluminium solid was the
    #         slave of the compliant 2 mm laminate).  Swapped, and the shell
    #         thickness is now accounted for so the bond actually closes.
    for nm, sl, ma in (('Int-1', 'WING-R-BOND', 'HINGE-A-BOND'),
                       ('Int-2', 'WING-L-BOND', 'HINGE-B-BOND')):
        m.SurfaceToSurfaceContactStd(
            name=nm, createStepName='Initial',
            master=SRF[ma], slave=SRF[sl], sliding=SMALL,
            thickness=ON, interactionProperty='IntProp-1',
            adjustMethod=OVERCLOSED, initialClearance=OMIT,
            datumAxis=None, clearanceRegion=None)
else:
    for nm, sl, ma in (('Bond-Right', 'WING-R-BOND', 'HINGE-A-BOND'),
                       ('Bond-Left',  'WING-L-BOND', 'HINGE-B-BOND')):
        m.Tie(name=nm, master=SRF[ma], slave=SRF[sl],
              positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON,
              thickness=ON)

# =========================== SECTION 8 - STEP ================================
# FIX #10: 100 increments (the default) could never complete a nlgeom step with
#          cohesive damage; limit raised and cut-back behaviour relaxed.
m.StaticStep(name='Step-1', previous='Initial', nlgeom=ON,
             timePeriod=1.0, maxNumInc=1000, initialInc=0.05,
             minInc=1e-08, maxInc=0.1,
             stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
             stabilizationMagnitude=2e-04, adaptiveDampingRatio=0.05,
             continueDampingFactors=True)
m.steps['Step-1'].control.setValues(
    allowPropagation=OFF, resetDefaultValues=OFF, discontinuous=ON,
    timeIncrementation=(8.0, 10.0, 9.0, 16.0, 10.0, 4.0, 12.0, 10.0, 6.0, 3.0, 50.0))

# FIX #11: the SMOOTH STEP amplitude existed in the deck but was never used by
#          any load - the loads ramped linearly.  It is applied here.
m.SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((0.0, 0.0), (1.0, 1.0)))

# ===================== SECTION 9 - BOUNDARY CONDITIONS =======================
# BC-1: pin nodes restrained in 2 and 3 (radial), free to rotate/slide axially
m.DisplacementBC(name='BC-1-Pin', createStepName='Initial',
                 region=a.sets['PIN-BC-NODES'], u1=UNSET, u2=SET, u3=SET,
                 ur1=UNSET, ur2=UNSET, ur3=UNSET, distributionType=UNIFORM)
# BC-2..BC-4: reference points - free in 2, 3 and rotation about 1 (hinge axis)
for nm, rp in (('BC-2-PinRP', 'RP-PIN-CENTER'),
               ('BC-3-RootRight', 'RP-HINGE-RIGHT'),
               ('BC-4-RootLeft', 'RP-HINGE-LEFT')):
    m.DisplacementBC(name=nm, createStepName='Initial', region=a.sets[rp],
                     u1=SET, u2=UNSET, u3=UNSET,
                     ur1=UNSET, ur2=SET, ur3=SET, distributionType=UNIFORM)

# =========================== SECTION 10 - LOADS ==============================
F = TIP_FORCE * LOAD_SCALE
M = ROOT_MOMENT * LOAD_SCALE
m.ConcentratedForce(name='Load-1-LeftWingUp',  createStepName='Step-1',
                    region=a.sets['RP-WING-LEFT'],  cf2=F, amplitude='Amp-1',
                    distributionType=UNIFORM, follower=ON)
m.ConcentratedForce(name='Load-2-LeftTipDown', createStepName='Step-1',
                    region=a.sets['RP-TIP-LEFT'],    cf2=-F, amplitude='Amp-1',
                    distributionType=UNIFORM, follower=ON)
m.ConcentratedForce(name='Load-3-RightWingUp', createStepName='Step-1',
                    region=a.sets['RP-WING-RIGHT'], cf2=F, amplitude='Amp-1',
                    distributionType=UNIFORM, follower=ON)
m.ConcentratedForce(name='Load-4-RightTipDown', createStepName='Step-1',
                    region=a.sets['RP-TIP-RIGHT'],   cf2=-F, amplitude='Amp-1',
                    distributionType=UNIFORM, follower=OFF)
m.Moment(name='Load-5-RootMomentRight', createStepName='Step-1',
         region=a.sets['RP-HINGE-RIGHT'], cm2=-M, amplitude='Amp-1',
         distributionType=UNIFORM)
m.Moment(name='Load-6-RootMomentLeft',  createStepName='Step-1',
         region=a.sets['RP-HINGE-LEFT'],  cm2=M, amplitude='Amp-1',
         distributionType=UNIFORM)

# ======================= SECTION 11 - OUTPUT REQUESTS ========================
# FIX #12: PE / PEEQ / PEMAG were requested although no material was plastic.
#          Plasticity now exists, so these variables are meaningful; damage
#          (DAMAGEFT...) and status (STATUS) are added for the failure study.
nodeVars = ('U', 'RF', 'CF')
elemVars = ['S', 'LE', 'PE', 'PEEQ', 'PEMAG',
            'HSNFTCRT', 'HSNFCCRT', 'HSNMTCRT', 'HSNMCCRT']
if PROGRESSIVE_DAMAGE:
    elemVars += ['DAMAGEFT', 'DAMAGEFC', 'DAMAGEMT', 'DAMAGEMC', 'STATUS', 'SDEG']
m.FieldOutputRequest(name='F-Output-1', createStepName='Step-1',
                     variables=tuple(nodeVars) + tuple(elemVars) +
                               ('CSTRESS', 'CDISP'),
                     numIntervals=50)
m.HistoryOutputRequest(name='H-Output-1', createStepName='Step-1',
                       variables=PRESELECT, numIntervals=200)

# ============================ SECTION 12 - JOB ===============================
if JOB_NAME in mdb.jobs.keys():
    del mdb.jobs[JOB_NAME]
job = mdb.Job(name=JOB_NAME, model=MODEL_NAME, type=ANALYSIS,
              description='l12 - hinged composite wing (revised)',
              memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
              nodalOutputPrecision=SINGLE,
              numCpus=N_CPUS, numGPUs=0,
              multiprocessingMode=DEFAULT, resultsFormat=ODB)
job.writeInput(consistencyChecking=OFF)
print('Input written: %s.inp' % JOB_NAME)

try:
    vp = session.viewports[session.currentViewportName]
    vp.setValues(displayedObject=a)
    vp.assemblyDisplay.setValues(mesh=ON, loads=ON, bcs=ON, constraints=ON,
                                 interactions=ON)
    vp.view.setValues(session.views['Iso']); vp.view.fitView()
except Exception as e:
    print('viewport skipped: %s' % e)

print('Model build complete.')
print('Input deck written: %s.inp' % JOB_NAME)
print('SUBMIT_JOB = %s' % SUBMIT_JOB)

if SUBMIT_JOB:
    job.submit(consistencyChecking=OFF)
    print('Job submitted on %d cpus' % N_CPUS)
    if WAIT_FOR_JOB:
        job.waitForCompletion()
        print('Job status: %s' % job.status)
        odb_file = os.path.join(os.getcwd(), JOB_NAME + '.odb')
        if os.path.isfile(odb_file):
            odb = session.openOdb(name=odb_file)
            try:
                vp = session.viewports[session.currentViewportName]
                vp.setValues(displayedObject=odb)
                vp.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
                vp.odbDisplay.setPrimaryVariable(variableLabel='S',
                        outputPosition=INTEGRATION_POINT,
                        refinement=(INVARIANT, 'Mises'))
                vp.view.fitView()
            except Exception as e:
                print('plot skipped: %s' % e)
            if 'Step-1' in odb.steps.keys() and len(odb.steps['Step-1'].frames):
                fr = odb.steps['Step-1'].frames[-1]
                umax = max((v.magnitude for v in fr.fieldOutputs['U'].values), default=0.0)
                smax = max((v.mises for v in fr.fieldOutputs['S'].values), default=0.0)
                print('t = %.3f | max |U| = %.3f mm | max Mises = %.1f MPa' %
                      (fr.frameValue, umax, smax))
            else:
                print('ODB opened, but Step-1 contains no frames.')
print('Done.')
