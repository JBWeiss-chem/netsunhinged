import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const $ = (id) => document.getElementById(id);
const S = {
  st: null, model: null, opts: null, colorMode: 'element', sc: [2, 2, 2],
  result: null, jobId: null, anim: null, frame: null, sel: null, regime: null,
};
// opened structures (see "loaded structures"): the one on screen and the one being read
const entries = [];
let active = null, pending = null;
const PALETTE = ['#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f', '#edc948', '#b07aa1',
  '#ff9da7', '#9c755f', '#bab0ac', '#86bcb6', '#d37295', '#fabfd2', '#8cd17d', '#499894', '#f1ce63'];
// result colours (Okabe & Ito, colour-blind safe): positive greens for flexible / semi-flexible,
// negative vermillion for rigid; semi-flexible lines are also dashed so colour is not needed
const CLS_COLOR = { flexible: '#009E73', 'semi-flexible': '#4FB892', rigid: '#D55E00' };
const CLS_DASH = { flexible: '', 'semi-flexible': '5 3', rigid: '' };

// ------------------------------------------------------------------ API
async function api(path, body) {
  const r = await fetch(path, body ? {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body),
  } : {});
  const j = await r.json();
  if (j.error) throw Object.assign(new Error(j.error), { status: r.status });
  return j;
}
function fail(e) { console.error(e); alert(e.message || e); }
// every text that comes from a structure file, the server or the catalogue is escaped before
// it is put into HTML (defence in depth; the server also sanitises CIF text)
const ESC = { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;', '`': '&#96;' };
const esc = (v) => String(v).replace(/[&<>"'`]/g, (c) => ESC[c]);

// ------------------------------------------------------------------ three.js
const viewer = $('viewer');
// no preserveDrawingBuffer: the PNG export renders and reads the image in one step
const renderer = new THREE.WebGLRenderer({ antialias: true });
// HiDPI screens: at most 2 device pixels per CSS pixel (4K laptops report up to 3); large scenes
// drop to 1 while the view is dragged and return to full sharpness on release
const MAX_PIXEL_RATIO = 2;
const fullPixelRatio = () => Math.min(window.devicePixelRatio || 1, MAX_PIXEL_RATIO);
renderer.setPixelRatio(fullPixelRatio());
viewer.appendChild(renderer.domElement);
// The scene is drawn only when something has changed (camera, content, colours, size, animation),
// not on every display refresh, so an idle viewer costs no GPU time.
let needsRender = true;
function requestRender() { needsRender = true; }
const scene = new THREE.Scene();
scene.background = new THREE.Color('#11151d');
// Two cameras: the perspective slider dolly-zooms the perspective camera (keeping the framed
// height constant) and switches to a true orthographic camera at its left end, where cell
// edges stack exactly when looking along a lattice vector.
const perspCam = new THREE.PerspectiveCamera(35, 1, 0.1, 50000);
const orthoCam = new THREE.OrthographicCamera(-10, 10, 10, -10, -50000, 50000);
let camera = perspCam;
camera.position.set(30, 25, 40);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = false;           // direct, unsmoothed response to the mouse
controls.addEventListener('change', requestRender);
scene.add(new THREE.AmbientLight(0xffffff, 0.55));
const keyLight = new THREE.DirectionalLight(0xffffff, 1.6);
camera.add(keyLight);
keyLight.position.set(5, 8, 10);
scene.add(perspCam, orthoCam);
const root = new THREE.Group();
scene.add(root);

function viewHeight() {
  if (camera.isOrthographicCamera) return (orthoCam.top - orthoCam.bottom) / orthoCam.zoom;
  const dist = camera.position.distanceTo(controls.target);
  return 2 * dist * Math.tan(THREE.MathUtils.degToRad(perspCam.fov / 2));
}
function aspect() { return viewer.clientWidth / Math.max(viewer.clientHeight, 1); }
function placeCamera(h, dir, up) {
  // dir: unit vector from target to camera; h: visible height at the target
  const p = +$('persp').value;
  const tgt = controls.target;
  const prev = camera;
  if (p <= 0.001) {
    camera = orthoCam;
    orthoCam.top = h / 2; orthoCam.bottom = -h / 2;
    orthoCam.left = -h * aspect() / 2; orthoCam.right = h * aspect() / 2;
    orthoCam.zoom = 1;
    camera.position.copy(tgt).addScaledVector(dir, Math.max(h, 50));
  } else {
    camera = perspCam;
    perspCam.fov = 1 + p * 59;
    perspCam.aspect = aspect();
    const dist = h / (2 * Math.tan(THREE.MathUtils.degToRad(perspCam.fov / 2)));
    perspCam.near = Math.max(0.05, dist / 200);
    perspCam.far = dist + 50 * h;
    camera.position.copy(tgt).addScaledVector(dir, dist);
  }
  camera.up.copy(up);
  camera.lookAt(tgt);
  camera.updateProjectionMatrix();
  if (prev !== camera) { camera.add(keyLight); controls.object = camera; }
  controls.update();
  requestRender();
}
$('persp').oninput = () => {
  const dir = camera.position.clone().sub(controls.target).normalize();
  placeCamera(viewHeight(), dir, camera.up.clone());
};

// ---------- crystallographic axis compass (bottom right)
// compass drawn in a single colour that contrasts with the background (white on dark, black on white)
const compassInk = () => ($('whiteBg').checked ? '#000000' : '#ffffff');
const AX_COLORS = ['#ffffff', '#ffffff', '#ffffff'];
const axScene = new THREE.Scene();
const axCam = new THREE.OrthographicCamera(-1.45, 1.45, 1.45, -1.45, -10, 10);
const axArrows = AX_COLORS.map((c) => {
  const arr = new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), new THREE.Vector3(), 1, c, 0.3, 0.17);
  arr.line.material.linewidth = 2;
  axScene.add(arr);
  return arr;
});
function labelTexture(text, color) {
  const cv = document.createElement('canvas'); cv.width = cv.height = 128;
  const g = cv.getContext('2d');
  g.font = 'italic bold 92px Georgia, "Times New Roman", serif'; g.textAlign = 'center'; g.textBaseline = 'middle';
  g.fillStyle = color; g.fillText(text, 64, 70);
  return new THREE.CanvasTexture(cv);
}
const axLabels = ['a', 'b', 'c'].map((t) => {
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: labelTexture(t, '#ffffff'), depthTest: false }));
  sp.scale.set(0.62, 0.62, 1);
  axScene.add(sp);
  return sp;
});
function setCompassInk() {
  const ink = compassInk();
  axArrows.forEach((a) => a.setColor(ink));
  axLabels.forEach((sp, i) => { sp.material.map.dispose(); sp.material.map = labelTexture('abc'[i], ink); sp.material.needsUpdate = true; });
  requestRender();
}
function updateCompass(L) {
  L.forEach((v, i) => {
    const d = new THREE.Vector3(...v).normalize();
    axArrows[i].setDirection(d);
    axLabels[i].position.copy(d).multiplyScalar(1.28);
  });
}
function renderCompass() {
  const size = 150, w = viewer.clientWidth;
  // stay clear of the animation / path bar at the bottom of the viewer
  const bar = [$('frameBar'), $('animBar')].find((b) => !b.hidden);
  const y0 = bar ? viewer.getBoundingClientRect().bottom - bar.getBoundingClientRect().top + 6 : 6;
  $('legend').style.bottom = `${bar ? y0 + 2 : 8}px`;
  axCam.position.copy(camera.position).sub(controls.target).normalize().multiplyScalar(5);
  axCam.up.copy(camera.up);
  axCam.lookAt(0, 0, 0);
  renderer.autoClear = false;
  renderer.clearDepth();
  renderer.setScissorTest(true);
  renderer.setScissor(w - size - 6, y0, size, size);
  renderer.setViewport(w - size - 6, y0, size, size);
  renderer.render(axScene, axCam);
  renderer.setScissorTest(false);
  renderer.setViewport(0, 0, w, viewer.clientHeight);
  renderer.autoClear = true;
}

const sphereGeo = new THREE.SphereGeometry(1, 20, 14);
const cylGeo = new THREE.CylinderGeometry(1, 1, 1, 10, 1);
// large scenes: atoms are then only a few pixels wide, so coarser shapes look the same and
// cut the triangle count about 4-fold (spheres) and 1.7-fold (bonds)
const LARGE_SCENE = 8000;                 // drawn atoms
const sphereGeoCoarse = new THREE.SphereGeometry(1, 10, 7);
const cylGeoCoarse = new THREE.CylinderGeometry(1, 1, 1, 6, 1);
let largeScene = false;
const material = new THREE.MeshStandardMaterial({ roughness: 0.45, metalness: 0.05 });
const jointMat = new THREE.MeshStandardMaterial({ color: '#CC79A7', roughness: 0.4, transparent: true, opacity: 0.9 });
const ballMat = new THREE.MeshStandardMaterial({ color: '#CC79A7', roughness: 0.4, transparent: true, opacity: 0.3, depthWrite: false });
let atomMesh = null, bondMesh = null, hingeMesh = null, ballMesh = null, cellLines = null;
let atomInst = [], bondInst = [], hingeInst = [], ballInst = [];
// topology view: nodes/linkers as spheres, net edges as two half-cylinders in unit colours.
// topoT blends between the atomistic (0) and the topology (1) picture, also during animations.
let topoEdgeMesh = null, topoEdgeInst = [];
let topoT = 0;

new ResizeObserver(() => {
  const w = viewer.clientWidth, h = viewer.clientHeight;
  renderer.setPixelRatio(fullPixelRatio());   // e.g. after moving the window to another screen
  renderer.setSize(w, h);
  requestRender();
  perspCam.aspect = w / Math.max(h, 1);
  perspCam.updateProjectionMatrix();
  const hh = orthoCam.top - orthoCam.bottom;
  orthoCam.left = -hh * perspCam.aspect / 2; orthoCam.right = hh * perspCam.aspect / 2;
  orthoCam.updateProjectionMatrix();
}).observe(viewer);

// small vector helpers on plain arrays
const vadd = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const vsub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const vscale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const imgVec = (n, L) => [0, 1, 2].map((k) => n[0] * L[0][k] + n[1] * L[1][k] + n[2] * L[2][k]);
const matmul = (A, B) => A.map((r) => [0, 1, 2].map((k) => r[0] * B[0][k] + r[1] * B[1][k] + r[2] * B[2][k]));

function atomVisible(s) {
  // dangling non-H groups (methyl C, terminal nitriles) are always drawn
  const a = S.st.atoms[s];
  if (a.role === 'guest') return $('showGuests').checked;
  if (a.el === 'H' || a.el === 'D') return $('showH').checked;
  return true;
}
function atomRadius(s) {
  const a = S.st.atoms[s];
  const r = (a.el === 'H' ? 0.17 : 0.2 + 0.22 * a.r) * (1 - topoT);
  // ticked atoms are drawn slightly larger so they stand out, not only by colour
  return S.model && $('showDof').checked && a.role === 'core' && S.model.flex[s] ? r * 1.35 : r;
}

function currentGeometry() {
  let X, L;
  if (S.frame) { X = S.frame.xyz; L = S.frame.L; } else { X = S.st.atoms.map((a) => a.xyz); L = S.st.L; }
  if (S.anim) {
    const f = S.anim.amp * Math.sin(S.anim.t * 2.4);
    X = X.map((x, s) => vadd(x, vscale(S.anim.disp[s], f)));
    if (S.anim.eps) {
      const LE = matmul(L, S.anim.eps);
      L = L.map((r, i) => vadd(r, vscale(LE[i], f)));
    }
  }
  return { X, L };
}

function images() {
  const out = [];
  for (let i = 0; i < S.sc[0]; i++) for (let j = 0; j < S.sc[1]; j++) for (let k = 0; k < S.sc[2]; k++) out.push([i, j, k]);
  return out;
}

function disposeMesh(m) { if (m) { root.remove(m); m.dispose?.(); } }

function buildScene() {
  [atomMesh, bondMesh, hingeMesh, ballMesh, topoEdgeMesh].forEach(disposeMesh);
  if (cellLines) { root.remove(cellLines); cellLines.geometry.dispose(); }
  atomMesh = bondMesh = hingeMesh = ballMesh = topoEdgeMesh = null;
  requestRender();
  if (!S.st) return;
  const imgs = images();
  atomInst = []; bondInst = []; hingeInst = []; ballInst = []; topoEdgeInst = [];
  const inside = (m) => m.every((v, k) => v >= 0 && v < S.sc[k]);
  const ghosts = new Set();
  const addGhost = (s, m) => {
    const key = `${s}:${m}`;
    if (!ghosts.has(key)) { ghosts.add(key); atomInst.push({ s, img: m }); }
  };
  // a bond segment runs from atom s (image ns) to the midpoint with atom t (image nt)
  const addBond = (k, s, ns, t, nt) => {
    bondInst.push({ b: k, s, ns, t, nt });
    bondInst.push({ b: k, s: t, ns: nt, t: s, nt: ns });
  };
  const boundary = $('boundaryAtoms').checked;
  for (const img of imgs) {
    S.st.atoms.forEach((a, s) => { if (atomVisible(s)) atomInst.push({ s, img }); });
    S.st.bonds.forEach((b, k) => {
      if (!(atomVisible(b.i) && atomVisible(b.j))) return;
      const nj = vadd(img, b.n), ni = vsub(img, b.n);
      if (inside(nj) || boundary) addBond(k, b.i, img, b.j, nj);
      else bondInst.push({ b: k, s: b.i, ns: img, t: b.j, nt: nj });
      if (!inside(nj) && boundary) addGhost(b.j, nj);
      if (!inside(ni)) {
        if (boundary) { addBond(k, b.j, img, b.i, ni); addGhost(b.i, ni); }
        else bondInst.push({ b: k, s: b.j, ns: img, t: b.i, nt: ni });
      }
    });
    buildTopology(img, inside, boundary);
    if (S.model && $('showJoints').checked && !S.anim && !S.frame) {
      S.model.joints.forEach((j) => {
        if (j.kind === 'hinge') hingeInst.push({ j, img }); else if (j.kind === 'ball') ballInst.push({ j, img });
      });
    }
  }
  largeScene = atomInst.length > LARGE_SCENE;
  atomMesh = new THREE.InstancedMesh(largeScene ? sphereGeoCoarse : sphereGeo, material, Math.max(atomInst.length, 1));
  bondMesh = new THREE.InstancedMesh(largeScene ? cylGeoCoarse : cylGeo, material, Math.max(bondInst.length, 1));
  atomMesh.count = atomInst.length; bondMesh.count = bondInst.length;
  root.add(atomMesh, bondMesh);
  if (hingeInst.length) { hingeMesh = new THREE.InstancedMesh(cylGeo, jointMat, hingeInst.length); root.add(hingeMesh); }
  if (ballInst.length) { ballMesh = new THREE.InstancedMesh(sphereGeo, ballMat, ballInst.length); root.add(ballMesh); }
  topoEdgeMesh = new THREE.InstancedMesh(cylGeo, material, Math.max(topoEdgeInst.length, 1));
  topoEdgeMesh.count = topoEdgeInst.length;
  root.add(topoEdgeMesh);
  cellLines = new THREE.LineSegments(new THREE.BufferGeometry(),
    new THREE.LineBasicMaterial({ color: $('whiteBg').checked ? '#444' : '#9aa5b8' }));
  cellLines.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(72), 3));
  root.add(cellLines);
  updateMatrices();
  updateColors();
}

const _m = new THREE.Matrix4(), _q = new THREE.Quaternion(), _p = new THREE.Vector3(), _s = new THREE.Vector3();
const _Y = new THREE.Vector3(0, 1, 0), _d = new THREE.Vector3();
function setCyl(mesh, idx, a, b, r) {
  _d.set(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
  const len = _d.length() || 1e-6;
  _q.setFromUnitVectors(_Y, _d.multiplyScalar(1 / len));
  _p.set((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2);
  _s.set(r, len, r);
  mesh.setMatrixAt(idx, _m.compose(_p, _q, _s));
}
function setSphere(mesh, idx, x, r) {
  _p.set(x[0], x[1], x[2]); _q.identity(); _s.set(r, r, r);
  mesh.setMatrixAt(idx, _m.compose(_p, _q, _s));
}

// ------------------------------------------------------------------ topology view
// lines only: every net edge is two half-cylinders, each in the colour of its building block
// (the block's representative atom, i.e. its most common non-H element, coloured by the
// current colour setting)
function buildTopology(img, inside, boundary) {
  const T = S.st.topology;
  // an edge segment runs from unit a (image na) to the midpoint with unit b (image nb)
  const both = (a, na, b, nb) => { topoEdgeInst.push({ a, na, b, nb }); topoEdgeInst.push({ a: b, na: nb, b: a, nb: na }); };
  T.edges.forEach(([u, v, t]) => {
    const nv = vadd(img, t), nu = vsub(img, t);
    if (inside(nv) || boundary) both(u, img, v, nv); else topoEdgeInst.push({ a: u, na: img, b: v, nb: nv });
    if (!inside(nu)) {
      if (boundary) both(v, img, u, nu); else topoEdgeInst.push({ a: v, na: img, b: u, nb: nu });
    }
  });
}
function unitPositions(X, L) {
  return S.st.topology.units.map((u) => {
    const p = [0, 0, 0];
    u.members.forEach(([s, m]) => { const x = vadd(X[s], imgVec(m, L)); p[0] += x[0]; p[1] += x[1]; p[2] += x[2]; });
    return vscale(p, 1 / u.members.length);
  });
}
function unitColor(u) { return baseAtomColor(S.st.topology.units[u].rep); }
function updateTopology(X, L) {
  if (!topoEdgeMesh) return;
  const P = unitPositions(X, L);
  topoEdgeInst.forEach((it, k) => {
    const a = vadd(P[it.a], imgVec(it.na, L)), b = vadd(P[it.b], imgVec(it.nb, L));
    setCyl(topoEdgeMesh, k, a, vscale(vadd(a, b), 0.5), 0.2 * topoT);
  });
  topoEdgeMesh.instanceMatrix.needsUpdate = true;
  topoEdgeMesh.computeBoundingSphere();
}

function updateMatrices() {
  if (!atomMesh) return;
  const { X, L } = currentGeometry();
  atomInst.forEach((it, k) => setSphere(atomMesh, k, vadd(X[it.s], imgVec(it.img, L)), atomRadius(it.s)));
  bondInst.forEach((it, k) => {
    const a = vadd(X[it.s], imgVec(it.ns, L));
    const c = vscale(vadd(a, vadd(X[it.t], imgVec(it.nt, L))), 0.5);
    setCyl(bondMesh, k, a, c, 0.11 * (1 - topoT));
  });
  atomMesh.instanceMatrix.needsUpdate = true;
  bondMesh.instanceMatrix.needsUpdate = true;
  if (hingeMesh) {
    hingeInst.forEach((it, k) => {
      const o = imgVec(it.img, L), [p0, p1] = it.j.points;
      const d = vsub(p1, p0), len = Math.hypot(...d) || 1, ext = vscale(d, 0.55 / len);
      setCyl(hingeMesh, k, vadd(vsub(p0, ext), o), vadd(vadd(p1, ext), o), 0.07 * (1 - topoT));
    });
    hingeMesh.instanceMatrix.needsUpdate = true;
  }
  if (ballMesh) {
    ballInst.forEach((it, k) => setSphere(ballMesh, k, vadd(it.j.points[0], imgVec(it.img, L)), 0.62 * (1 - topoT)));
    ballMesh.instanceMatrix.needsUpdate = true;
  }
  updateTopology(X, L);
  atomMesh.visible = bondMesh.visible = topoT < 1;
  topoEdgeMesh.visible = topoT > 0;
  // unit cell box
  const pos = cellLines.geometry.attributes.position.array;
  const corners = [];
  for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) for (let k = 0; k < 2; k++) corners.push(imgVec([i, j, k], L));
  const edges = [[0, 1], [2, 3], [4, 5], [6, 7], [0, 2], [1, 3], [4, 6], [5, 7], [0, 4], [1, 5], [2, 6], [3, 7]];
  edges.forEach(([a, b], e) => { pos.set(corners[a], e * 6); pos.set(corners[b], e * 6 + 3); });
  cellLines.geometry.attributes.position.needsUpdate = true;
  cellLines.geometry.computeBoundingSphere();
  updateCompass(L);
  atomMesh.computeBoundingSphere(); bondMesh.computeBoundingSphere();
  requestRender();
}

// ------------------------------------------------------------------ colours
function isFlex(s) {
  const a = S.st.atoms[s];
  return S.model ? S.model.flex[s] : (S.opts && S.opts.flex_angle_labels.includes(a.label));
}
function bondState(k) {
  if (!S.model) return S.st.bonds[k].auto;
  return S.model.bond_state[String(k)] || null;
}
function bondKey(k) {
  const b = S.st.bonds[k];
  return [S.st.atoms[b.i].label, S.st.atoms[b.j].label].sort().join('|');
}
// colour-blind safe (Okabe & Ito) highlight colours
const HL = { flex: '#E69F00', hinge: '#56B4E9', sel: '#F0E442', joint: '#CC79A7' };
function dofOn() { return $('showDof').checked && S.model; }
function baseAtomColor(s) {
  const a = S.st.atoms[s];
  if (S.colorMode === 'role') return a.role === 'core' ? '#748ffc' : a.role === 'pruned' ? '#adb5bd' : '#f06595';
  if (S.colorMode === 'body') {
    if (!S.model || S.model.site_body[s] < 0) return '#555';
    return PALETTE[S.model.site_body[s] % PALETTE.length];
  }
  return a.color;
}
function atomColor(s) {
  const a = S.st.atoms[s];
  if (S.sel && S.sel.type === 'atom' && a.label === S.st.atoms[S.sel.s].label && a.role === 'core') return HL.sel;
  if (dofOn() && a.role === 'core' && isFlex(s)) return HL.flex;
  return baseAtomColor(s);
}
function bondColor(k, s) {
  const b = S.st.bonds[k];
  if (S.sel && S.sel.type === 'bond' && b.core && bondKey(k) === S.sel.key) return HL.sel;
  if (dofOn() && b.core && bondState(k) === 'hinge' && S.model.bond_status[bondKey(k)] !== 'sbu') return HL.hinge;
  return baseAtomColor(s);
}
const _c = new THREE.Color();
function updateColors() {
  if (!atomMesh) return;
  atomInst.forEach((it, k) => atomMesh.setColorAt(k, _c.set(atomColor(it.s))));
  bondInst.forEach((it, k) => bondMesh.setColorAt(k, _c.set(bondColor(it.b, it.s))));
  if (atomMesh.instanceColor) atomMesh.instanceColor.needsUpdate = true;
  if (bondMesh.instanceColor) bondMesh.instanceColor.needsUpdate = true;
  topoEdgeInst.forEach((it, k) => topoEdgeMesh.setColorAt(k, _c.set(unitColor(it.a))));
  if (topoEdgeMesh.instanceColor) topoEdgeMesh.instanceColor.needsUpdate = true;
  renderLegend();
  requestRender();
}
function renderLegend() {
  const L = $('legend');
  const item = (c, t) => `<i style="background:${c}"></i>${t}`;
  if (!S.st) { L.innerHTML = ''; return; }
  let html;
  if ($('topology').checked) {
    const seen = new Map();
    S.st.topology.units.forEach((u, k) => { if (!seen.has(u.label)) seen.set(u.label, [u, k]); });
    L.innerHTML = 'topology: ' + [...seen.values()].map(([u, k]) =>
      item(unitColor(k), `${u.kind === 'node' ? 'node' : 'linker'} ${esc(u.label)}`)).join('');
    return;
  }
  if (S.colorMode === 'role') html = item('#748ffc', 'framework (analysed)') +
    item('#adb5bd', 'dangling (H, terminal)') + item('#f06595', 'guest / solvent');
  else if (S.colorMode === 'body') html = `${S.model ? S.model.bodies.length : 0} rigid bodies`;
  else html = [...new Set(S.st.atoms.map((a) => a.el))].map((e) =>
    item(S.st.atoms.find((a) => a.el === e).color, esc(e))).join('');
  if (dofOn()) html += '<br>' + item(HL.flex, 'flexible angles') + item(HL.hinge, 'rotatable bond');
  if ($('showJoints').checked && S.model && !S.anim && !S.frame) html += item(HL.joint, 'hinge axis / ball joint');
  L.innerHTML = html;
}

// ------------------------------------------------------------------ camera
function frameCamera(dir) {
  if (!S.st) return;
  const L = S.st.L;
  const sc = S.sc;
  const center = vscale(vadd(vadd(imgVec([sc[0], 0, 0], L), imgVec([0, sc[1], 0], L)), imgVec([0, 0, sc[2]], L)), 0.5);
  const diag = Math.hypot(...imgVec(sc, L));
  const d = dir ? dir : [0.6, 0.45, 1];
  const n = Math.hypot(...d);
  controls.target.set(...center);
  const up = new THREE.Vector3(0, 1, 0);
  if (dir && Math.abs(dir[1]) / n > 0.95) up.set(0, 0, 1);
  placeCamera(1.15 * diag, new THREE.Vector3(...d).normalize(), up);
}

// ------------------------------------------------------------------ picking
const ray = new THREE.Raycaster();
const mouse = new THREE.Vector2();
function pick(ev) {
  if (!atomMesh) return null;
  const r = renderer.domElement.getBoundingClientRect();
  mouse.set(((ev.clientX - r.left) / r.width) * 2 - 1, -((ev.clientY - r.top) / r.height) * 2 + 1);
  ray.setFromCamera(mouse, camera);
  const hits = ray.intersectObjects([atomMesh, bondMesh], false);
  if (!hits.length) return null;
  const h = hits[0];
  if (h.object === atomMesh) return { type: 'atom', s: atomInst[h.instanceId].s };
  const bi = bondInst[h.instanceId].b;
  return { type: 'bond', b: bi, key: bondKey(bi) };
}
let downPos = null;
renderer.domElement.addEventListener('pointerdown', (e) => { downPos = [e.clientX, e.clientY]; });
renderer.domElement.addEventListener('pointerup', (e) => {
  if (!downPos || Math.hypot(e.clientX - downPos[0], e.clientY - downPos[1]) > 4) return;
  const p = pick(e);
  S.sel = p;
  renderSelection();
  updateColors();
});
let hoverT = 0;
let dragLowRes = false;
renderer.domElement.addEventListener('pointermove', (e) => {
  if (e.buttons) {
    // rotating / panning: no hover picking (it tests every atom and bond), and large scenes
    // render at 1 device pixel per CSS pixel until the button is released
    $('tooltip').hidden = true;
    if (largeScene && !dragLowRes && fullPixelRatio() > 1) {
      dragLowRes = true; renderer.setPixelRatio(1); requestRender();
    }
    return;
  }
  const now = performance.now();
  if (now - hoverT < 60) return;
  hoverT = now;
  const p = pick(e);
  const tt = $('tooltip');
  if (!p) { tt.hidden = true; return; }
  const r = viewer.getBoundingClientRect();
  tt.style.left = `${e.clientX - r.left + 12}px`; tt.style.top = `${e.clientY - r.top + 12}px`;
  if (p.type === 'atom') {
    const a = S.st.atoms[p.s];
    tt.textContent = `${a.label}  ${a.env}${isFlex(p.s) ? '  · flexible angles' : ''}`;
  } else {
    const b = S.st.bonds[p.b];
    tt.textContent = `${p.key.replace('|', '–')}  ${b.len.toFixed(3)} Å  ${b.core ? bondState(p.b) : 'not analysed'}`;
  }
  tt.hidden = false;
});
renderer.domElement.addEventListener('pointerleave', () => { $('tooltip').hidden = true; });
// the button may be released outside the viewer
window.addEventListener('pointerup', () => {
  if (dragLowRes) { dragLowRes = false; renderer.setPixelRatio(fullPixelRatio()); requestRender(); }
});

// ------------------------------------------------------------------ selection panel
function labelsWithEnv(env) {
  return Object.entries(S.st.labels).filter(([, d]) => d.env === env).map(([l]) => l);
}
function renderSelection() {
  const box = $('selection');
  if (!S.sel || !S.opts) { box.hidden = true; return; }
  box.hidden = false;
  let html = '<span class="close" id="selClose">✕</span>';
  if (S.sel.type === 'atom') {
    const s = S.sel.s, a = S.st.atoms[s], info = S.st.labels[a.label];
    html += `<h4>${esc(a.label)} <span class="muted">${esc(a.env)}</span></h4>`;
    html += `<div class="muted">${a.role === 'core' ? `framework atom · ${info.count} symmetry copies` :
      a.role === 'pruned' ? 'dangling atom (removed for analysis)' : 'guest / solvent (removed)'}
      ${a.sbu ? '· SBU' : ''}${a.bridging ? ' · carboxylate-type donor' : ''}
      ${S.model && S.model.site_body[s] >= 0 ? `· body ${S.model.site_body[s]}` : ''}</div>`;
    if (a.role === 'core') {
      const labOn = S.opts.flex_angle_labels.includes(a.label);
      const envLabels = labelsWithEnv(a.env);
      const envOn = envLabels.every((l) => S.opts.flex_angle_labels.includes(l));
      const siteOn = S.opts.flex_angle_sites.includes(s);
      html += `<div style="margin-top:6px"><b>Bond angles</b> at this atom:
        ${isFlex(s) ? '<span style="color:#b07400">flexible (pivot)</span>' : 'rigid'}</div>
        <div class="btns">
          <button class="small ${labOn ? 'on' : ''}" data-act="lab">label ${esc(a.label)} (${info.count})</button>
          <button class="small ${envOn ? 'on' : ''}" data-act="env">all ${esc(a.env)} (${envLabels.length} labels)</button>
          <button class="small ${siteOn ? 'on' : ''}" data-act="site">this atom only</button>
        </div>`;
    }
  } else {
    const k = S.sel.b, b = S.st.bonds[k];
    const [l1, l2] = S.sel.key.split('|');
    const ov = S.opts.bond_overrides[S.sel.key];
    html += `<h4>${esc(l1)}–${esc(l2)} <span class="muted">${b.len.toFixed(3)} Å</span></h4>`;
    if (!b.core) html += '<div class="muted">involves a removed atom – not part of the analysis</div>';
    else {
      html += `<div class="muted">automatic: ${esc(b.auto)} (${esc(b.reason)})</div>
        <div style="margin-top:6px"><b>Torsion</b>: ${bondState(k)}</div>
        <div class="btns">
          <button class="small ${!ov ? 'on' : ''}" data-act="auto">auto</button>
          <button class="small ${ov === 'hinge' ? 'on' : ''}" data-act="hinge">hinge (rotatable)</button>
          <button class="small ${ov === 'locked' ? 'on' : ''}" data-act="locked">locked</button>
        </div>`;
    }
  }
  box.innerHTML = html;
  $('selClose').onclick = () => { S.sel = null; renderSelection(); updateColors(); };
  box.querySelectorAll('button[data-act]').forEach((btn) => { btn.onclick = () => selectionAction(btn.dataset.act); });
}
function toggleIn(arr, v) { const i = arr.indexOf(v); if (i >= 0) arr.splice(i, 1); else arr.push(v); }
function selectionAction(act) {
  const o = S.opts;
  if (S.sel.type === 'atom') {
    const a = S.st.atoms[S.sel.s];
    if (act === 'lab') toggleIn(o.flex_angle_labels, a.label);
    if (act === 'site') toggleIn(o.flex_angle_sites, S.sel.s);
    if (act === 'env') {
      const labs = labelsWithEnv(a.env);
      const allOn = labs.every((l) => o.flex_angle_labels.includes(l));
      labs.forEach((l) => { const i = o.flex_angle_labels.indexOf(l); if (allOn && i >= 0) o.flex_angle_labels.splice(i, 1); if (!allOn && i < 0) o.flex_angle_labels.push(l); });
    }
  } else {
    if (act === 'auto') delete o.bond_overrides[S.sel.key]; else o.bond_overrides[S.sel.key] = act;
  }
  optionsChanged();
}

// ------------------------------------------------------------------ DOF tables
function renderTables() {
  const aStat = S.model ? S.model.atom_status : {}, bStat = S.model ? S.model.bond_status : {};
  const SBU = '<span class="pill sbu" title="fixed by the rigid SBU">SBU</span>';
  // metal-related rows first, then alphabetical by element
  const atomRel = (d) => d.relevant;
  const lt = $('labelTable').querySelector('tbody');
  const labs = Object.entries(S.st.labels).sort((x, y) =>
    (atomRel(y[1]) - atomRel(x[1])) || (x[1].el + x[0]).localeCompare(y[1].el + y[0], undefined, { numeric: true }));
  lt.innerHTML = labs.map(([l, d]) => {
    const st = aStat[l], on = S.opts.flex_angle_labels.includes(l);
    // non-metal SBU atoms (e.g. mu-O) are fully governed by the rigid SBU; metals can release it
    const locked = st === 'sbu' && !d.metal;
    const cell = locked ? SBU : `<input type="checkbox" ${on ? 'checked' : ''}>${st === 'sbu' ? ' ' + SBU : ''}`;
    return `<tr data-label="${esc(l)}" class="${atomRel(d) ? '' : 'minor'} ${on ? 'on-flex' : ''}"><td>${esc(l)}</td>
      <td>${esc(d.env)}${d.bridging ? ' <span class="muted" title="carboxylate-type donor">⚲</span>' : ''}</td>
      <td>${d.count}</td><td class="nowrap">${cell}</td></tr>`;
  }).join('');
  lt.querySelectorAll('input').forEach((inp) => {
    inp.onchange = () => { toggleIn(S.opts.flex_angle_labels, inp.closest('tr').dataset.label); optionsChanged(); };
  });
  lt.querySelectorAll('tr').forEach((tr) => {
    tr.onmouseenter = () => { const s = S.st.atoms.findIndex((a) => a.label === tr.dataset.label && a.role === 'core'); if (s >= 0) { S.sel = { type: 'atom', s }; updateColors(); } };
    tr.onmouseleave = () => { S.sel = null; updateColors(); renderSelection(); };
  });
  const bt = $('bondTable').querySelector('tbody');
  const bcs = Object.entries(S.st.bond_classes).sort((x, y) =>
    (y[1].relevant - x[1].relevant) || x[1].elements.localeCompare(y[1].elements) || x[0].localeCompare(y[0], undefined, { numeric: true }));
  bt.innerHTML = bcs.map(([k, d]) => {
    const ov = S.opts.bond_overrides[k];
    const st = bStat[k] || (ov || d.auto);
    const rot = st === 'hinge';
    const cell = st === 'sbu' ? SBU : `<input type="checkbox" ${rot ? 'checked' : ''}>` +
      (ov ? ' <a class="reset" title="back to automatic (' + esc(d.auto) + ')">↺</a>' : '');
    return `<tr data-key="${esc(k)}" class="${d.relevant ? '' : 'minor'} ${rot && st !== 'sbu' ? 'on-hinge' : ''}">
      <td>${esc(k.replace('|', '–'))}</td><td>${d.count}</td>
      <td title="${esc(d.reason)}">${d.auto === 'hinge' ? 'rot.' : d.auto}</td><td class="nowrap">${cell}</td></tr>`;
  }).join('');
  bt.querySelectorAll('input').forEach((inp) => {
    inp.onchange = () => {
      const k = inp.closest('tr').dataset.key, want = inp.checked ? 'hinge' : 'locked';
      // store an override only when it differs from the automatic choice
      if (want === S.st.bond_classes[k].auto) delete S.opts.bond_overrides[k]; else S.opts.bond_overrides[k] = want;
      optionsChanged();
    };
  });
  bt.querySelectorAll('a.reset').forEach((a) => {
    a.onclick = () => { delete S.opts.bond_overrides[a.closest('tr').dataset.key]; optionsChanged(); };
  });
  bt.querySelectorAll('tr').forEach((tr) => {
    tr.onmouseenter = () => { S.sel = { type: 'bond', key: tr.dataset.key, b: S.st.bonds.findIndex((b, k) => b.core && bondKey(k) === tr.dataset.key) }; updateColors(); };
    tr.onmouseleave = () => { S.sel = null; updateColors(); renderSelection(); };
  });
  // the checkbox shows the effective state: metals ticked by hand release their SBU
  const metalLabels = Object.entries(S.st.labels).filter(([, d]) => d.metal).map(([l]) => l);
  const nFlexMetals = metalLabels.filter((l) => S.opts.flex_angle_labels.includes(l)).length;
  const box = $('rigidSbu');
  box.checked = S.opts.rigid_sbu && nFlexMetals === 0;
  box.indeterminate = S.opts.rigid_sbu && nFlexMetals > 0 && nFlexMetals < metalLabels.length;
  $('sbuCutoff').value = S.opts.sbu_cutoff;
  updateControlStates();
}

let modelTimer = null;
function optionsChanged() {
  S.modelVersion = (S.modelVersion || 0) + 1;
  S.edited = true;
  S.modelDirty = true;
  renderPresetInfo();
  renderTables();
  renderSelection();
  clearTimeout(modelTimer);
  modelTimer = setTimeout(() => refreshModel().catch(fail), 200);
}

function renderPresetInfo() {
  const p = S.preset;
  if (!p) { $('presetInfo').innerHTML = ''; return; }
  $('presetInfo').innerHTML = S.edited
    ? `Current model: <b>custom</b> (edited from “${esc(p.label)}”)`
    : `Current model: <b>${esc(p.label)}</b>${p.reason ? `<br><span class="muted">${esc(p.reason)}</span>` : ''}`;
}
async function refreshModel(presetName) {
  if (!S.st) return;
  if (presetName) S.modelVersion = (S.modelVersion || 0) + 1;
  const st = S.st, version = S.modelVersion;
  $('modelInfo').textContent = 'building model…';
  const model = await withSession(() => api('/api/model', presetName
    ? { id: S.st.id, preset: presetName } : { id: S.st.id, options: S.opts }));
  // the user switched to another structure meanwhile: it is rebuilt when shown again
  if (S.st !== st) return;
  if (S.modelVersion === version) S.modelDirty = false;
  S.model = model;
  S.opts = S.model.options;
  if (S.model.preset) { S.preset = S.model.preset; S.edited = false; }
  renderPresetInfo();
  invalidateResults();
  renderTables();
  renderModelInfo();
  renderInfInfo();
  updateControlStates();
  renderSelection();
  if (!S.anim && !S.frame) buildScene(); else updateColors();
}

function invalidateResults() {
  S.result = null;
  clearResults();
  renderEntries();
}
function clearResults() {
  S.regime = null;
  $('verdict').innerHTML = ''; $('plot').innerHTML = '';
  $('regTable').querySelector('tbody').innerHTML = '';
  $('exportRes').hidden = true;
  $('plotView').hidden = true;
  stopFrame(false);
}

// ------------------------------------------------------------------ info panels
function renderStructInfo() {
  const st = S.st;
  const c = st.cell;
  $('structInfo').innerHTML = `<div class="kv">
    <span>name</span><span>${esc(st.name)}</span>
    <span>formula <i class="info" tabindex="0" data-tip="Framework formula without guests (solvent, counter-ions and other atoms not bonded to the framework), reduced to the smallest formula unit; Z is the number of formula units per unit cell. H atoms are included as far as they are present in the CIF."></i></span><span>${esc(st.framework_formula_reduced)}${st.z > 1 ? ` <span class="muted">(Z = ${st.z})</span>` : ''}</span>
    <span>cell</span><span>${c.slice(0, 3).map((x) => x.toFixed(3)).join(', ')} Å<br>${c.slice(3).map((x) => x.toFixed(2)).join(', ')}°</span>
    <span>atoms</span><span>${st.n_atoms} (analysed ${st.n_core}, dangling ${st.n_pruned}, guests ${st.n_guest})</span>
    <span>SBUs / cell</span><span>${st.n_sbu}</span>
    <span>periodicity</span><span>${st.dimensionality}D bonded network</span></div>`;
  $('warnings').innerHTML = st.warnings.map((w) => `<div class="warn">${esc(w)}</div>`).join('');
  const rg = $('removeGuests');
  rg.hidden = !st.n_guest;
  rg.textContent = `Remove guests (${st.n_guest} atoms: ${guestFormula()})`;
}
function guestFormula() {
  const c = {};
  S.st.atoms.forEach((a) => { if (a.role === 'guest') c[a.el] = (c[a.el] || 0) + 1; });
  return Object.entries(c).map(([e, n]) => `${e}${n > 1 ? n : ''}`).join(' ');
}
$('removeGuests').onclick = async () => {
  try {
    stopFrame(false); S.anim = null; $('animBar').hidden = true; S.sel = null; renderSelection();
    // framework atoms keep their order, so per-site choices can be remapped
    const newIndex = {};
    let k = 0;
    S.st.atoms.forEach((a, s) => { if (a.role !== 'guest') newIndex[s] = k++; });
    S.opts.flex_angle_sites = S.opts.flex_angle_sites.filter((s) => s in newIndex).map((s) => newIndex[s]);
    const st = S.st;
    const reduced = await withSession(() => api('/api/remove_guests', { id: S.st.id }));
    if (S.st !== st) return;
    S.st = reduced; S.guestsRemoved = true;
    renderStructInfo();
    // keep the user's degrees of freedom: atom labels are unchanged
    await refreshModel();
  } catch (e) { fail(e); }
};
function renderModelInfo() {
  const s = S.model.summary;
  $('modelInfo').innerHTML = `<div class="kv">
    <span>rigid bodies</span><span>${s.n_bodies}${s.n_linear_bodies ? ` (${s.n_linear_bodies} linear)` : ''}</span>
    <span>hinges</span><span>${s.n_hinges}</span>
    <span>ball joints</span><span>${s.n_ball_joints}</span>
    <span>pivot atoms</span><span>${s.n_flex_angle_atoms}</span>
    <span>bonds hinge / locked</span><span>${s.n_hinge_bonds} / ${s.n_locked_bonds}</span>
    <span>joint constraints</span><span>${s.n_constraints}</span></div>
    ${s.n_sbu_released ? `<p class="hint">${s.n_sbu_released} SBU(s) are not rigid because their metals have flexible angles.</p>` : ''}`;
}
function renderInfInfo() {
  const r = S.model.rigidity;
  if (r.skipped) {
    $('infInfo').innerHTML = `<p class="hint">This model has ${r.n_dof} unknowns, more than this server's limit for the infinitesimal analysis. The finite deformation test below still works.</p>`;
    return;
  }
  let html = `<div class="kv"><span>mechanisms</span><span>${r.n_mechanisms}</span>
    <span>cell-changing</span><span><b>${r.n_macroscopic}</b></span>
    <span>internal (local)</span><span>${r.n_internal}</span></div>`;
  if (r.n_macroscopic === 0) html += '<p class="hint">The model allows no deformation of the unit cell without distorting a rigid body or opening a joint: it is infinitesimally rigid for this choice of degrees of freedom.</p>';
  else html += '<p class="hint">Cell deformations the model allows to first order (no rigid body distorted, no joint opened). In the material these are soft directions, not free ones. Follow them to finite amplitude with the deformation test.</p>';
  r.macroscopic_modes.forEach((m, k) => {
    html += `<div class="mode"><span class="k">Mode ${k + 1}</span> · ${esc(m.kind)}<br>
      <span class="muted">${esc(m.description)} · ΔV/V ${m.volume_rel >= 0 ? '+' : ''}${m.volume_rel.toFixed(2)} · strain share ${(100 * m.strain_share).toFixed(m.strain_share < 0.1 ? 1 : 0)}%</span>
      <div><button class="small" data-mode="${k}">▶ animate</button></div></div>`;
  });
  if (r.internal_modes.length) {
    html += `<div style="margin-top:6px">Internal mechanisms: ${r.internal_modes.map((_, k) =>
      `<button class="small" data-imode="${k}">${k + 1}</button>`).join(' ')}</div>`;
  }
  // internal modes that still strain the cell a little (bent linkers in experimental structures)
  if (r.n_internal_weak_strain) html += `<p class="hint">${r.n_internal_weak_strain} internal mechanism(s) change the cell slightly, but less than ${100 * r.strain_share_min}% of their motion comes from the cell strain: they are local motions (typically ring spins of slightly bent linkers) and are not followed as cell-changing modes.</p>`;
  $('infInfo').innerHTML = html;
  $('infInfo').querySelectorAll('[data-mode]').forEach((b) => {
    b.onclick = () => {
      const m = r.macroscopic_modes[+b.dataset.mode], mx = m.max_disp || 1;
      startAnim(`Mode ${+b.dataset.mode + 1}: ${m.kind}`, m.disp.map((d) => vscale(d, 1 / mx)),
        m.eps.map((row) => vscale(row, 1 / mx)));
    };
  });
  $('infInfo').querySelectorAll('[data-imode]').forEach((b) => {
    b.onclick = () => {
      const m = r.internal_modes[+b.dataset.imode];
      startAnim(`Internal mechanism ${+b.dataset.imode + 1}`, m.disp, m.strain_share > 1e-6 ? m.eps : null);
    };
  });
}

// ------------------------------------------------------------------ animation
function startAnim(title, disp, eps) {
  S.frame = null; $('frameBar').hidden = true;
  S.anim = { disp, eps, t: 0, amp: +$('animAmp').value };
  $('animTitle').textContent = title;
  $('animBar').hidden = false;
  buildScene();
  updateControlStates();
}
function stopAnim() { S.anim = null; $('animBar').hidden = true; buildScene(); updateControlStates(); }
$('animStop').onclick = stopAnim;
$('animAmp').oninput = () => { if (S.anim) S.anim.amp = +$('animAmp').value; };

let lastT = performance.now();
function loop(t) {
  const dt = (t - lastT) / 1000; lastT = t;
  // smooth cross-fade between atomistic and topology view
  const target = $('topology').checked ? 1 : 0;
  const blending = topoT !== target;
  if (blending) topoT = target > topoT ? Math.min(1, topoT + dt / 0.4) : Math.max(0, topoT - dt / 0.4);
  if (S.anim) S.anim.t += dt;
  if (S.st && (S.anim || blending)) updateMatrices();        // requests a render
  if (needsRender) {
    needsRender = false;
    renderer.render(scene, camera);
    if (S.st) renderCompass();
  }
  requestAnimationFrame(loop);
}
requestAnimationFrame(loop);
// safety net: any use of a control may change the picture (and the bars the compass avoids)
['input', 'change', 'click'].forEach((t) => document.addEventListener(t, requestRender, true));

// ------------------------------------------------------------------ finite deformation
$('run').onclick = async () => {
  try {
    invalidateResults();
    const settings = {
      dmin: +$('dmin').value, dmax: +$('dmax').value, step: +$('step').value, follow_modes: $('followModes').checked,
      extended: $('extended').checked, paper_regimes: true,
    };
    // the job belongs to this structure: it keeps running when another structure is shown,
    // and its results are stored with this one
    const entry = active;
    const state = () => (entry === active ? S : entry.saved);
    const { job } = await withSession(() => api('/api/analyze', { id: S.st.id, options: S.opts, settings }));
    const own = state();
    if (!own) return;
    own.jobId = job;
    own.activeJob = job;         // only the latest job of a structure may touch its state
    own.progress = 0;
    // results belong to this structure and model; edits during the run make them stale
    const stamp = { st: own.st.id, version: own.modelVersion };
    showProgress(); renderEntries();
    const poll = async () => {
      if (state()?.activeJob !== job) return;       // superseded, abandoned or closed
      const j = await api(`/api/job?job=${job}`);
      const cur = state();
      if (cur?.activeJob !== job) return;
      cur.progress = j.progress;
      if (j.status === 'queued' || j.status === 'running') {
        if (entry === active) {
          showProgress();
          $('verdict').innerHTML = j.status === 'queued' ? '<p class="hint">Waiting for a free worker (other analyses are running)…</p>' : '';
        }
        setTimeout(poll, j.status === 'queued' ? 1500 : 400);
        return;
      }
      cur.activeJob = null;
      if (entry === active) showProgress();
      renderEntries();
      if (j.status === 'error') { fail(new Error(entry === active ? j.message : `${entry.label}: ${j.message}`)); return; }
      if (!cur.st || cur.st.id !== stamp.st || cur.modelVersion !== stamp.version) {
        if (entry === active) $('verdict').innerHTML = '<p class="hint">The structure or the DOFs were changed during the run, so its results were discarded. Run the analysis again.</p>';
        return;
      }
      cur.result = j.result; cur.result.settings = settings;
      if (entry === active) renderResults();
      renderEntries();
    };
    poll();
  } catch (e) { showProgress(); fail(e); }
};
// the run button and progress bar show the job of the structure on screen
function showProgress() {
  const running = !!(S.st && S.activeJob);
  $('run').disabled = !S.st || running;
  $('progress').hidden = !running;
  $('progress').firstElementChild.style.width = `${Math.round((S.progress || 0) * 100)}%`;
}

// The plotted metric and the thresholds are view settings: the server returns E and F for every
// point, and the regimes are classified here from F, so both can be changed after the run.
function metricKey() { return $('metric').value === 'absolute' ? 'E' : 'F'; }
function thresholds() {
  const v = (id, def) => { const x = +$(id).value; return Number.isFinite(x) && x > 0 ? x : def; };
  return [v('thr', 0.01), v('thr2', 0.05)];
}
const clsColor = (c) => CLS_COLOR[c];

// same criteria as deform.classify (keep the two in sync); F is scale-free, absolute E is not,
// so the classification always uses F, whichever metric is plotted
function classify(F, d, thr, thr2) {
  if (F.every((v) => v < thr)) return ['flexible', 'mismatch below threshold over the whole range'];
  const below = (side) => {
    const ks = d.map((x, i) => i).filter((i) => side * (d[i] - 1) > 1e-9);
    return ks.length > 0 && ks.every((i) => F[i] < thr);
  };
  if (below(-1)) return ['semi-flexible', 'mismatch below threshold for d < 1 only (one-sided)'];
  if (below(1)) return ['semi-flexible', 'mismatch below threshold for d > 1 only (one-sided)'];
  const i1 = d.reduce((b, x, i) => (Math.abs(x - 1) < Math.abs(d[b] - 1) ? i : b), 0);
  let lo = i1, hi = i1;
  while (lo - 1 >= 0 && F[lo - 1] < thr) lo--;
  while (hi + 1 < F.length && F[hi + 1] < thr) hi++;
  if (d[hi] - d[lo] >= 0.039) return ['semi-flexible', `mismatch below threshold only for ${d[lo].toFixed(2)} <= d <= ${d[hi].toFixed(2)}`];
  if (F.every((v) => v < thr2)) return ['semi-flexible', 'low but non-zero mismatch over the whole range'];
  return ['rigid', ''];
}
function classifyResult(R) {
  const [thr, thr2] = thresholds();
  Object.assign(R.settings, { thr, thr2 });
  Object.values(R.regimes).forEach((r) => { [r.classification, r.detail] = classify(r.F, r.d, thr, thr2); });
  const cls = Object.values(R.regimes).map((r) => r.classification);
  R.verdict = cls.includes('flexible') ? 'flexible' : cls.includes('semi-flexible') ? 'semi-flexible' : 'rigid';
  return R;
}

function renderResults() {
  const R = classifyResult(S.result), key = metricKey();
  const counts = { flexible: 0, 'semi-flexible': 0, rigid: 0 };
  Object.values(R.regimes).forEach((r) => counts[r.classification]++);
  const flexList = Object.entries(R.regimes).filter(([, r]) => r.classification === 'flexible').map(([k]) => k);
  $('verdict').innerHTML = `<span class="big" style="color:${CLS_COLOR[R.verdict]}">${R.verdict === 'flexible' ? 'Potentially flexible' : R.verdict === 'semi-flexible' ? 'Semi-flexible' : 'Rigid'}</span>
    <div class="muted">${counts.flexible} flexible · ${counts['semi-flexible']} semi · ${counts.rigid} rigid regimes · ${R.seconds}s
    ${flexList.length ? `<br>flexible: ${esc(flexList.join(', '))}` : ''}</div>`;
  $('maxHead').textContent = `max ${key}`;
  const tb = $('regTable').querySelector('tbody');
  tb.innerHTML = Object.entries(R.regimes).map(([k, r]) => {
    const mx = Math.max(...r[key]);
    // volume range spanned by the scan, relative to the original cell
    const V0 = cellVolume(S.st.cell), rel = r.cell.map((c) => 100 * (cellVolume(c) / V0 - 1));
    const fmt = (v) => `${v >= 0 ? '+' : ''}${v.toFixed(Math.abs(v) >= 10 ? 0 : 1)}`;
    return `<tr class="click${k === S.regime ? ' active' : ''}" data-reg="${esc(k)}" title="${esc(r.description)}${r.detail ? ' — ' + esc(r.detail) : ''}">
      <td><b>${esc(k)}</b> <span class="muted">${esc(r.description)}</span>${r.detail ? `<div class="det">${esc(r.detail)}</div>` : ''}</td>
      <td><span class="tag ${r.classification}">${r.classification === 'semi-flexible' ? 'semi' : r.classification}</span></td>
      <td>${mx.toExponential(1)}</td><td>${fmt(Math.min(...rel))}…${fmt(Math.max(...rel))}</td></tr>`;
  }).join('');
  tb.querySelectorAll('tr').forEach((tr) => { tr.onclick = () => selectRegime(tr.dataset.reg).catch(fail); });
  $('exportRes').hidden = false;
  $('plotView').hidden = false;
  renderPlot();
}
// a new metric or threshold only redraws: nothing is recomputed
['metric', 'thr', 'thr2'].forEach((id) => {
  $(id).addEventListener('input', () => {
    if (!S.result) return;
    renderResults();
    if (S.frame) { frameInfo(S.frame); renderGraphs(); }
    renderEntries();
  });
});

function renderPlot() {
  const R = S.result, key = metricKey(), st = R.settings;
  const W = 330, H = 210, ml = 38, mr = 8, mt = 8, mb = 26;
  const lo = -10, hi = 1;
  const x = (d) => ml + ((d - st.dmin) / (st.dmax - st.dmin)) * (W - ml - mr);
  const y = (v) => { const l = Math.min(hi, Math.max(lo, Math.log10(Math.max(v, 1e-30)))); return mt + ((hi - l) / (hi - lo)) * (H - mt - mb); };
  let svg = `<svg viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg">`;
  for (let e = lo; e <= hi; e += 2) svg += `<line x1="${ml}" x2="${W - mr}" y1="${y(10 ** e)}" y2="${y(10 ** e)}" stroke="#eee"/><text x="${ml - 4}" y="${y(10 ** e) + 3}" text-anchor="end">${e === lo ? '≤' : ''}1e${e}</text>`;
  for (let d = Math.ceil(st.dmin * 20 - 1e-9) / 20; d <= st.dmax + 1e-9; d += 0.05) svg += `<text x="${x(d)}" y="${H - 12}" text-anchor="middle">${d.toFixed(2)}</text>`;
  svg += `<text x="${(W + ml) / 2}" y="${H - 1}" text-anchor="middle">order parameter d</text>`;
  svg += `<text x="10" y="${mt + 4}" >${key === 'F' ? 'F' : 'E'}</text>`;
  svg += `<line x1="${x(1)}" x2="${x(1)}" y1="${mt}" y2="${H - mb}" stroke="#ccc"/>`;
  if (key === 'F') {
    [[st.thr, CLS_COLOR.flexible], [st.thr2, CLS_COLOR['semi-flexible']]].forEach(([t, c]) => {
      svg += `<line x1="${ml}" x2="${W - mr}" y1="${y(t)}" y2="${y(t)}" stroke="${c}" stroke-dasharray="1.5 2.5" stroke-width="1.4"/>`;
    });
  }
  // E(1) = 0 exactly (drawn at the bottom of the axis). F(1) = 0/0 is undefined, but F is continuous
  // there (E and E_affine both grow as (d - 1)^2), so the two branches are joined across d = 1.
  const curve = (r) => r.d.map((d, i) => [d, r[key][i]]).filter(([d]) => key === 'E' || Math.abs(d - 1) > 1e-9);
  const pts = (c) => c.map(([d, v]) => `${x(d)},${y(v)}`).join(' ');
  const entries = Object.entries(R.regimes);
  entries.forEach(([k, r]) => {
    if (k === S.regime) return;
    const col = S.regime ? '#c3c9d2' : clsColor(r.classification);
    svg += `<polyline fill="none" stroke="${col}" stroke-dasharray="${S.regime ? '' : CLS_DASH[r.classification] || ''}" stroke-width="${S.regime ? 0.8 : 1.1}" opacity="${S.regime ? 0.7 : 0.9}"
      points="${pts(curve(r))}"><title>${esc(k)}</title></polyline>`;
  });
  // the selected regime on top: halo, thick line, data points, the frame shown and a label
  const r = S.regime && R.regimes[S.regime];
  if (r) {
    const c = curve(r), col = clsColor(r.classification);
    svg += `<polyline fill="none" stroke="#fff" stroke-width="5" stroke-linejoin="round" points="${pts(c)}"/>`;
    svg += `<polyline fill="none" stroke="${col}" stroke-width="2.4" stroke-linejoin="round" points="${pts(c)}"><title>${esc(S.regime)}</title></polyline>`;
    c.forEach(([d, v]) => { svg += `<circle cx="${x(d)}" cy="${y(v)}" r="2.3" fill="${col}" stroke="#fff" stroke-width="0.8"/>`; });
    if (S.frame && Math.abs(S.frame.d - 1) > 1e-9) {
      svg += `<circle cx="${x(S.frame.d)}" cy="${y(S.frame[key])}" r="5" fill="none" stroke="#222" stroke-width="1.4"/>`;
    }
    svg += `<text x="${W - mr - 2}" y="${mt + 10}" text-anchor="end" class="sel" fill="${col}">${esc(S.regime)} · ${r.classification}</text>`;
  }
  svg += '</svg>';
  $('plot').innerHTML = svg;
}

async function selectRegime(name) {
  S.regime = name;
  document.querySelectorAll('#regTable tr').forEach((tr) => tr.classList.toggle('active', tr.dataset.reg === name));
  renderPlot();
  const r = S.result.regimes[name];
  const sl = $('frameSlider');
  sl.max = r.d.length - 1;
  sl.value = 0;
  $('frameTitle').textContent = name;
  $('frameBar').hidden = false;
  await loadFrame();
}

// ------------------------------------------------------------------ path graphs
const SERIES = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#000000', '#999999'];
function lineChart(title, unit, x, series, mark) {
  const W = 330, H = 150, ml = 40, mr = 8, mt = 18, mb = 24;
  if (!series.length) return `<div class="chart"><div class="ct">${title}</div><p class="hint">nothing changes along this path</p></div>`;
  const all = series.flatMap((s) => s.y);
  let lo = Math.min(0, ...all), hi = Math.max(0, ...all);
  if (hi - lo < 1e-6) { hi += 1; lo -= 1; }
  const pad = 0.06 * (hi - lo); lo -= pad; hi += pad;
  const x0 = Math.min(...x), x1 = Math.max(...x);
  const X = (v) => ml + ((v - x0) / (x1 - x0)) * (W - ml - mr);
  const Y = (v) => mt + ((hi - v) / (hi - lo)) * (H - mt - mb);
  const nice = (v) => (Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2));
  let svg = `<svg viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg">`;
  svg += `<text x="${ml}" y="11" class="ct">${title} <tspan class="u">[${unit}]</tspan></text>`;
  for (let k = 0; k <= 4; k++) {
    const v = lo + (k / 4) * (hi - lo);
    svg += `<line x1="${ml}" x2="${W - mr}" y1="${Y(v)}" y2="${Y(v)}" stroke="#eee"/><text x="${ml - 4}" y="${Y(v) + 3}" text-anchor="end">${nice(v)}</text>`;
  }
  svg += `<line x1="${ml}" x2="${W - mr}" y1="${Y(0)}" y2="${Y(0)}" stroke="#bbb"/>`;
  for (let d = Math.ceil(x0 * 20) / 20; d <= x1 + 1e-9; d += 0.05) svg += `<text x="${X(d)}" y="${H - 10}" text-anchor="middle">${d.toFixed(2)}</text>`;
  svg += `<text x="${(W + ml) / 2}" y="${H - 1}" text-anchor="middle">d</text>`;
  if (mark != null) svg += `<line x1="${X(mark)}" x2="${X(mark)}" y1="${mt}" y2="${H - mb}" stroke="#888" stroke-dasharray="3 3"/>`;
  series.forEach((s, i) => {
    const c = SERIES[i % SERIES.length];
    svg += `<polyline fill="none" stroke="${c}" stroke-width="1.6" points="${x.map((v, k) => `${X(v)},${Y(s.y[k])}`).join(' ')}"><title>${esc(s.name)}</title></polyline>`;
    x.forEach((v, k) => { svg += `<circle cx="${X(v)}" cy="${Y(s.y[k])}" r="1.6" fill="${c}"/>`; });
  });
  svg += '</svg>';
  const legend = series.map((s, i) => `<span><i style="background:${SERIES[i % SERIES.length]}"></i>${esc(s.name)}</span>`).join('');
  return `<div class="chart">${svg}<div class="cl">${legend}</div></div>`;
}
function renderGraphs() {
  if (!S.result || !S.regime) { $('graphCard').hidden = true; return; }
  const r = S.result.regimes[S.regime], g = r.geometry;
  const mark = S.frame ? S.frame.d : null;
  $('graphCard').hidden = false;
  $('graphTitle').textContent = S.regime;
  $('graphHint').textContent = `${r.description}. Changes relative to the original structure (d = 1). ` +
    'For each angle/torsion type the instance that changes most is shown (n = number of instances).';
  const key = metricKey();
  const html = [
    lineChart('Cell lengths and volume', '%', g.d, [
      { name: 'a', y: g.lengths.a }, { name: 'b', y: g.lengths.b }, { name: 'c', y: g.lengths.c },
      { name: 'V', y: g.volume }], mark),
    lineChart('Cell angles', 'Δ°', g.d, [
      { name: 'α', y: g.angles_cell.alpha }, { name: 'β', y: g.angles_cell.beta }, { name: 'γ', y: g.angles_cell.gamma }], mark),
    lineChart(`Flexible bond angles (${g.n_flex_angles})`, 'Δ°', g.d,
      g.flex_angles.map((c) => ({ name: `${c.name} (n=${c.n})`, y: c.y })), mark),
    lineChart(`Torsions about hinge bonds (${g.n_torsions})`, 'Δ°', g.d,
      g.torsions.map((c) => ({ name: `${c.name} (n=${c.n})`, y: c.y })), mark),
    lineChart(key === 'F' ? 'Relative mismatch F = E/E_affine' : 'Joint mismatch E', key === 'F' ? '–' : 'Å²/joint', g.d,
      [{ name: key, y: r[key] }], mark),
  ];
  $('regGraphs').innerHTML = html.join('');
}

// ------------------------------------------------------------------ trajectory export
const FORMATS = [['extxyz', 'extended XYZ (.xyz)'], ['pdb', 'multi-model PDB'], ['traj', 'ASE .traj']];
const FMT_EXT = { extxyz: 'xyz', pdb: 'pdb', traj: 'traj' };
['animFmt', 'pathFmt'].forEach((id) => {
  $(id).innerHTML = FORMATS.map(([v, t]) => `<option value="${v}">${t}</option>`).join('');
});
function saveExport(res, base, fmt) {
  if (res.base64) {
    const bin = Uint8Array.from(atob(res.base64), (c) => c.charCodeAt(0));
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([bin], { type: 'application/octet-stream' }));
    a.download = `${base}.${FMT_EXT[fmt]}`; a.click();
  } else download(`${base}.${FMT_EXT[fmt]}`, res.text, 'text/plain');
}
$('animExport').onclick = async () => {
  if (!S.anim) return;
  try {
    // one period of the oscillation, sampled like the on-screen animation
    const nf = 24, saved = S.anim.t, frames = [];
    for (let k = 0; k < nf; k++) {
      S.anim.t = (k / nf) * (2 * Math.PI) / 2.4;
      const { X, L } = currentGeometry();
      frames.push({ xyz: X, L, note: `frame=${k} amplitude=${S.anim.amp}A` });
    }
    S.anim.t = saved;
    const fmt = $('animFmt').value;
    const res = await withSession(() => api('/api/export', { id: S.st.id, frames, fmt }));
    saveExport(res, `${S.st.stem}_${$('animTitle').textContent.split(':')[0].replace(/\s+/g, '_')}`, fmt);
  } catch (e) { fail(e); }
};
$('pathExport').onclick = async () => {
  try {
    const fmt = $('pathFmt').value;
    const res = await api('/api/export', { id: S.st.id, job: S.jobId, regime: S.regime, fmt });
    saveExport(res, `${S.st.stem}_${S.regime.replace(/\s+/g, '_')}_path`, fmt);
  } catch (e) { fail(e); }
};
// Slider moves fire many input events. Only one /api/frame request runs at a time; moves during
// it just mark that the newest position must be fetched next, so responses can never arrive out
// of order. Frames already fetched are cached (the key includes the job id).
const frameCache = new Map();
const FRAME_CACHE_MAX = 200;
let frameBusy = false, frameAgain = false;
async function loadFrame() {
  if (frameBusy) { frameAgain = true; return; }
  frameBusy = true;
  try {
    do {
      frameAgain = false;
      if (!S.result || !S.regime) break;
      const job = S.jobId, regime = S.regime;
      const d = S.result.regimes[regime].d[+$('frameSlider').value];
      const key = `${job}|${regime}|${d}`;
      let f = frameCache.get(key);
      if (!f) {
        f = await api(`/api/frame?job=${job}&regime=${encodeURIComponent(regime)}&d=${d}`);
        if (frameCache.size >= FRAME_CACHE_MAX) frameCache.delete(frameCache.keys().next().value);
        frameCache.set(key, f);
      }
      // a newer slider position, another regime or a new run supersedes this frame
      if (!frameAgain && S.jobId === job && S.regime === regime) showFrame(f);
    } while (frameAgain);
  } finally {
    frameBusy = false;
  }
}
function showFrame(f) {
  if (S.anim) { S.anim = null; $('animBar').hidden = true; }
  const first = !S.frame;
  S.frame = f;
  frameInfo(f);
  if (first) buildScene(); else updateMatrices();
  renderPlot();
  renderGraphs();
  updateControlStates();
}
function frameInfo(f) {
  const V = cellVolume(f.cell), V0 = cellVolume(S.st.cell);
  const dV = (100 * (V / V0 - 1));
  $('frameInfo').innerHTML = `d = ${f.d.toFixed(3)} · ${metricKey()} = ${(metricKey() === 'F' ? f.F : f.E).toExponential(1)} · ` +
    `a,b,c = ${f.cell.slice(0, 3).map((v) => v.toFixed(2)).join(', ')} Å · α,β,γ = ${f.cell.slice(3).map((v) => v.toFixed(1)).join(', ')}° · ` +
    `<b>V = ${V.toFixed(1)} Å³ (${dV >= 0 ? '+' : ''}${dV.toFixed(1)} %)</b>`;
}
$('frameSlider').oninput = () => loadFrame().catch(fail);
function stopFrame(rebuild = true) {
  const had = !!S.frame;
  S.frame = null; $('frameBar').hidden = true;
  document.querySelectorAll('#regTable tr').forEach((tr) => tr.classList.remove('active'));
  $('graphCard').hidden = true;
  if (had && rebuild) buildScene();
  updateControlStates();
}
$('frameStop').onclick = () => { S.regime = null; stopFrame(); if (S.result) renderPlot(); };

// ------------------------------------------------------------------ export
function download(name, text, type = 'application/json') {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([text], { type }));
  a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
$('exportOpts').onclick = () => download(`${S.st?.stem || 'netsunhinged'}_dof.json`, JSON.stringify(S.opts, null, 2));
$('importOpts').onchange = async (e) => {
  const f = e.target.files[0]; if (!f) return;
  try {
    const o = JSON.parse(await f.text());
    S.opts = { flex_angle_labels: o.flex_angle_labels || [], flex_angle_sites: o.flex_angle_sites || [],
      bond_overrides: o.bond_overrides || {}, rigid_sbu: !!o.rigid_sbu, sbu_cutoff: o.sbu_cutoff || 5 };
    optionsChanged();
  } catch (err) { fail(err); }
  e.target.value = '';
};
$('exportRes').onclick = () => download(`${S.st.stem}_flexibility.json`, JSON.stringify({
  structure: S.st.name, cell: S.st.cell, options: S.opts, model: S.model.summary,
  infinitesimal: { n_macroscopic: S.model.rigidity.n_macroscopic, n_internal: S.model.rigidity.n_internal,
    modes: S.model.rigidity.macroscopic_modes.map((m) => ({ kind: m.kind, description: m.description, strain: m.strain6 })) },
  finite: S.result,
}, null, 1));
function cellVolume(p) {
  const [a, b, c] = p, [al, be, ga] = p.slice(3).map((x) => (x * Math.PI) / 180);
  return a * b * c * Math.sqrt(1 - Math.cos(al) ** 2 - Math.cos(be) ** 2 - Math.cos(ga) ** 2
    + 2 * Math.cos(al) * Math.cos(be) * Math.cos(ga));
}
function cellParams(L) {
  const n = (v) => Math.hypot(...v), dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const [a, b, c] = L.map(n);
  const ang = (u, v, lu, lv) => (Math.acos(dot(u, v) / (lu * lv)) * 180) / Math.PI;
  return [a, b, c, ang(L[1], L[2], b, c), ang(L[0], L[2], a, c), ang(L[0], L[1], a, b)];
}
function inv3(M) {
  const [a, b, c] = M[0], [d, e, f] = M[1], [g, h, i] = M[2];
  const A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g;
  const det = a * A + b * B + c * C;
  return [[A / det, -(b * i - c * h) / det, (b * f - c * e) / det],
    [B / det, (a * i - c * g) / det, -(a * f - c * d) / det],
    [C / det, -(a * h - b * g) / det, (a * e - b * d) / det]];
}
$('frameCif').onclick = () => {
  if (!S.frame) return;
  const L = S.frame.L, Li = inv3(L), p = cellParams(L);
  const lines = [`data_${S.st.name}_${S.regime}_d${S.frame.d.toFixed(3)}`.replace(/[^\w.]/g, '_'),
    ...['a', 'b', 'c'].map((k, i) => `_cell_length_${k} ${p[i].toFixed(5)}`),
    ...['alpha', 'beta', 'gamma'].map((k, i) => `_cell_angle_${k} ${p[i + 3].toFixed(4)}`),
    "_symmetry_space_group_name_H-M 'P 1'", 'loop_', '_symmetry_equiv_pos_as_xyz', 'x,y,z',
    'loop_', '_atom_site_label', '_atom_site_type_symbol', '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z'];
  S.frame.xyz.forEach((x, s) => {
    const f = [0, 1, 2].map((k) => { const v = x[0] * Li[0][k] + x[1] * Li[1][k] + x[2] * Li[2][k]; return v - Math.floor(v); });
    lines.push(`${S.st.atoms[s].label} ${S.st.atoms[s].el} ${f.map((v) => v.toFixed(6)).join(' ')}`);
  });
  download(`${S.st.stem}_${S.regime}_d${S.frame.d.toFixed(3)}.cif`, lines.join('\n') + '\n', 'chemical/x-cif');
};
$('savePng').onclick = () => {
  renderer.render(scene, camera);
  renderCompass();
  const a = document.createElement('a');
  a.href = renderer.domElement.toDataURL('image/png');
  a.download = `${S.st?.stem || 'structure'}.png`; a.click();
};

// ------------------------------------------------------------------ controls
$('colorMode').onchange = () => { S.colorMode = $('colorMode').value; updateColors(); };
// ---------- supercell: a x b x c repeats, limited by the number of drawn atoms
const MAX_DRAWN_ATOMS = 40000;
function visibleAtomCount() { return S.st ? S.st.atoms.filter((_, s) => atomVisible(s)).length : 0; }
function showSupercell(note) {
  ['scA', 'scB', 'scC'].forEach((id, k) => { $(id).value = S.sc[k]; });
  const n = visibleAtomCount() * S.sc[0] * S.sc[1] * S.sc[2];
  $('scInfo').innerHTML = note || (S.st ? `${n.toLocaleString()} atoms` : '');
}
function setSupercell(sc) {
  sc = sc.map((v) => Math.min(8, Math.max(1, Math.round(v) || 1)));
  const n = visibleAtomCount() * sc[0] * sc[1] * sc[2];
  if (n > MAX_DRAWN_ATOMS) {
    showSupercell(`<span class="err">${sc.join('×')} = ${n.toLocaleString()} atoms (max ${MAX_DRAWN_ATOMS.toLocaleString()})</span>`);
    return false;
  }
  const changed = sc.some((v, k) => v !== S.sc[k]);
  S.sc = sc;
  showSupercell();
  if (changed) { buildScene(); frameCamera(); }
  return true;
}
['scA', 'scB', 'scC'].forEach((id) => {
  $(id).onchange = () => setSupercell(['scA', 'scB', 'scC'].map((x) => +$(x).value));
});
$('scCube').onclick = () => {
  if (!S.st) return;
  const len = S.st.L.map((v) => Math.hypot(...v)), target = Math.max(...len);
  let sc = len.map((l) => Math.max(1, Math.round(target / l)));
  // shrink uniformly until the atom limit is respected
  while (!setSupercell(sc) && sc.some((v) => v > 1)) {
    const k = sc.indexOf(Math.max(...sc));
    sc = sc.map((v, i) => (i === k ? v - 1 : v));
  }
};
$('showJoints').onchange = buildScene;
$('showH').onchange = () => { buildScene(); showSupercell(); };
$('showGuests').onchange = () => { buildScene(); showSupercell(); };
$('boundaryAtoms').onchange = buildScene;
$('whiteBg').onchange = () => { scene.background.set($('whiteBg').checked ? '#ffffff' : '#11151d'); setCompassInk(); buildScene(); };
$('topology').onchange = () => { renderLegend(); updateControlStates(); };
document.querySelectorAll('[data-view]').forEach((b) => {
  b.onclick = () => { if (!S.st) return; const k = 'abc'.indexOf(b.dataset.view); frameCamera(S.st.L[k]); };
});
$('resetDefaults').onclick = () => refreshModel('paper').catch(fail);
$('resetZeolite').onclick = () => refreshModel('zeolite').catch(fail);
$('presetRigid').onclick = () => refreshModel('rigid').catch(fail);
$('showDof').onchange = buildScene;
$('rigidSbu').onchange = () => {
  // ticked: every cluster is one rigid unit; unticked: the angles inside the clusters
  // (at the metals and at atoms bonded only to metals, e.g. mu-O) become flexible
  const rigid = $('rigidSbu').checked;
  S.opts.rigid_sbu = rigid;
  Object.entries(S.st.labels).forEach(([l, d]) => {
    if (!d.sbu_internal) return;
    const i = S.opts.flex_angle_labels.indexOf(l);
    if (rigid && i >= 0) S.opts.flex_angle_labels.splice(i, 1);
    if (!rigid && i < 0) S.opts.flex_angle_labels.push(l);
  });
  optionsChanged();
};
$('sbuCutoff').onchange = () => { S.opts.sbu_cutoff = +$('sbuCutoff').value; optionsChanged(); };

// ------------------------------------------------------------------ loaded structures
// Every structure of every opened CIF (a file may hold several data blocks) gets a box in the
// Structure card. A structure is built on the server when it is first shown; afterwards its
// model, DOF choices and results are kept in the browser, so switching back to it is instant.
// The state of the structure on screen lives in S, that of the others in entry.saved.
const ENTRY_KEYS = ['st', 'model', 'opts', 'preset', 'edited', 'modelVersion', 'modelDirty',
  'result', 'jobId', 'activeJob', 'progress', 'guestsRemoved'];
function freshState(st) {
  return { st, model: null, opts: null, preset: null, edited: false, modelVersion: 0, modelDirty: true,
    result: null, jobId: null, activeJob: null, progress: 0, guestsRemoved: false };
}

// adds a box for each structure in the CIF and returns the first; a file opened again is not
// added twice
async function openCif(src) {
  const same = (e) => (src.example ? e.src.example === src.example
    : e.src.text === src.text && e.src.name === src.name);
  const known = entries.find(same);
  if (known) return known;
  const { name, blocks } = await api('/api/blocks', src);
  const file = name || src.name || 'structure';
  const added = blocks.map((b, block) => ({
    src, block, file, multi: blocks.length > 1, label: blocks.length > 1 ? b : file, saved: null,
  }));
  entries.push(...added);
  return added[0];
}

let switching = 0;
async function activate(entry) {
  if (entry === active) return;
  const token = ++switching;
  pending = entry; renderEntries();
  try {
    if (!entry.saved) entry.saved = freshState(await api('/api/load', { ...entry.src, block: entry.block }));
  } catch (e) {
    if (pending === entry) pending = null;
    closeEntry(entry);
    fail(e);
    return;
  }
  // another structure was chosen, or this one closed, while it was read
  if (token !== switching || !entries.includes(entry)) return;
  pending = null;
  leaveActive();
  if (active) active.saved = Object.fromEntries(ENTRY_KEYS.map((k) => [k, S[k]]));
  active = entry;
  Object.assign(S, entry.saved);
  await showActive();
}
// stops what belongs to the structure on screen; an edit not yet sent to the server stays
// marked (modelDirty) and is applied when the structure is shown again
function leaveActive() {
  stopFrame(false); S.anim = null; $('animBar').hidden = true; S.sel = null; renderSelection();
  clearTimeout(modelTimer);
}
async function showActive() {
  // keep the current repeats if they fit, otherwise the default 2x2x2, otherwise 1x1x1
  const fits = (sc) => visibleAtomCount() * sc[0] * sc[1] * sc[2] <= MAX_DRAWN_ATOMS;
  if (!fits(S.sc)) S.sc = fits([2, 2, 2]) ? [2, 2, 2] : [1, 1, 1];
  showSupercell();
  renderStructInfo();
  renderEntries();
  clearResults();
  showProgress();
  try {
    if (S.model && !S.modelDirty) {
      renderPresetInfo(); renderTables(); renderModelInfo(); renderInfInfo(); renderSelection();
      buildScene();
      if (S.result) renderResults();
    } else {
      clearModelPanels();
      buildScene();
      await refreshModel(S.model ? undefined : 'auto');
    }
  } catch (e) { fail(e); }
  frameCamera();
  updateControlStates();
}
function clearModelPanels() {
  $('presetInfo').innerHTML = '';
  $('modelInfo').textContent = '–'; $('infInfo').textContent = '–';
  ['labelTable', 'bondTable'].forEach((id) => { $(id).querySelector('tbody').innerHTML = ''; });
}
function showEmpty() {
  leaveActive();
  active = null;
  Object.assign(S, freshState(null));
  clearResults(); clearModelPanels(); buildScene();
  $('structInfo').textContent = 'Load a CIF file to start.';
  $('warnings').innerHTML = ''; $('removeGuests').hidden = true;
  showSupercell(); showProgress(); updateControlStates(); renderEntries();
}
function closeEntry(entry) {
  const k = entries.indexOf(entry);
  if (k < 0) return;
  entries.splice(k, 1);
  entry.saved = null;          // its running analysis, if any, is no longer followed
  if (entry === active) {
    const next = entries[k] || entries[k - 1];
    showEmpty();
    if (next) activate(next);
  }
  renderEntries();
}

function renderEntries() {
  const box = $('structList');
  box.hidden = !entries.length;
  box.replaceChildren(...entries.map((e) => {
    const st = e === active ? S : e.saved;
    const row = document.createElement('div');
    row.className = `sitem${e === active ? ' active' : ''}${e === pending ? ' loading' : ''}`;
    const pick = document.createElement('button');
    pick.className = 'pick';
    pick.title = e.multi ? `${e.label} (${e.file})` : e.label;
    const mark = document.createElement('span');
    mark.className = 'mark';
    if (st?.activeJob) { mark.classList.add('running'); mark.title = 'deformation test running'; }
    else if (st?.result) {
      const v = classifyResult(st.result).verdict;
      mark.style.background = clsColor(v);
      mark.title = `finite test: ${v}`;
    }
    const nm = document.createElement('span');
    nm.className = 'nm'; nm.textContent = e.label;
    pick.append(mark, nm);
    if (e.multi) {
      const f = document.createElement('span');
      f.className = 'file'; f.textContent = e.file;
      pick.append(f);
    }
    pick.onclick = () => activate(e);
    const x = document.createElement('button');
    x.className = 'x'; x.textContent = '×'; x.title = 'close'; x.setAttribute('aria-label', `close ${e.label}`);
    x.onclick = () => closeEntry(e);
    row.append(pick, x);
    return row;
  }));
  box.querySelector('.active')?.scrollIntoView({ block: 'nearest' });
}

// a structure's server session ends after an hour without use; it is then rebuilt from its CIF
// (which the browser keeps) and the call is repeated. DOF choices and results stay as they are.
async function withSession(call) {
  try { return await call(); } catch (e) {
    if (e.status !== 410 || !active || !S.st) throw e;
    const st = S.st;
    let fresh = await api('/api/load', { ...active.src, block: active.block });
    if (S.guestsRemoved) fresh = await api('/api/remove_guests', { id: fresh.id });
    if (S.st !== st) throw e;
    st.id = fresh.id;          // same atoms in the same order: only the session changes
    return call();
  }
}

async function openAndShow(sources) {
  let first = null;
  const errors = [];
  for (const [label, src] of sources) {
    try { first ||= await openCif(await src()); } catch (e) { errors.push(`${label}: ${e.message}`); }
  }
  renderEntries();
  if (errors.length) fail(new Error(errors.join('\n')));
  if (first) await activate(first);
}
$('file').onchange = async (e) => {
  const files = [...e.target.files];
  e.target.value = '';
  await openAndShow(files.map((f) => [f.name, async () => ({ text: await f.text(), name: f.name.replace(/\.cif$/i, '') })]));
};
$('examples').onchange = () => {
  const v = $('examples').value;
  $('examples').value = '';
  if (v) openAndShow([[v, async () => ({ example: v })]]);
};
api('/api/examples').then((r) => {
  $('examples').innerHTML += r.groups.map((g) => `<optgroup label="${esc(g.category)}">${
    g.items.map((e) => `<option value="${esc(e.file)}">${esc(e.label)}</option>`).join('')}</optgroup>`).join('');
}).catch(() => {});

// ------------------------------------------------------------------ info icons
// a small "i" next to headings/controls; its data-tip text appears in a floating box on hover.
// Each line of the tip (&#10; in the attribute) is a paragraph; a leading "Ticked:" or "Unticked:"
// is highlighted. Built from text nodes only, so a tip can never inject markup.
(() => {
  const tip = document.createElement('div');
  tip.id = 'infotip'; tip.hidden = true; document.body.appendChild(tip);
  const render = (text) => {
    tip.replaceChildren(...text.split('\n').map((line) => {
      const p = document.createElement('div');
      const m = line.match(/^(Ticked|Unticked):\s*/);
      if (m) {
        const b = document.createElement('b');
        b.className = 'tick'; b.textContent = `${m[1]}:`;
        p.append(b, ` ${line.slice(m[0].length)}`);
      } else p.textContent = line;
      return p;
    }));
  };
  const show = (el) => {
    render(el.dataset.tip); tip.hidden = false;
    const r = el.getBoundingClientRect(), w = tip.offsetWidth, h = tip.offsetHeight;
    let x = Math.min(window.innerWidth - w - 8, Math.max(8, r.left + r.width / 2 - w / 2));
    let y = r.bottom + 6;
    if (y + h > window.innerHeight - 8) y = r.top - h - 6;
    tip.style.left = `${x}px`; tip.style.top = `${y}px`;
  };
  document.querySelectorAll('i.info').forEach((el) => { el.tabIndex = 0; });
  document.addEventListener('mouseover', (e) => { const el = e.target.closest?.('i.info'); if (el) show(el); });
  document.addEventListener('mouseout', (e) => { if (e.target.closest?.('i.info')) tip.hidden = true; });
  document.addEventListener('focusin', (e) => { if (e.target.matches?.('i.info')) show(e.target); });
  document.addEventListener('focusout', (e) => { if (e.target.matches?.('i.info')) tip.hidden = true; });
  // icons sit inside <label>s: a click must not toggle the associated checkbox
  document.addEventListener('click', (e) => { if (e.target.closest?.('i.info')) e.preventDefault(); }, true);
})();
setCompassInk();

// ------------------------------------------------------------------ control states
// controls that cannot have any effect in the current state are greyed out and locked
function setEnabled(id, on) {
  const el = $(id);
  if (!el) return;
  el.disabled = !on;
  // dim the label only when every control in it is locked; a label shared with a still
  // usable control must stay normal
  const lab = el.closest('label');
  if (lab) {
    const ctrls = [...lab.querySelectorAll('input, select, button')];
    lab.classList.toggle('off', ctrls.length > 0 && ctrls.every((c) => c.disabled));
  }
}
function updateControlStates() {
  const loaded = !!S.st, topo = $('topology').checked, moving = !!(S.anim || S.frame);
  // nothing to act on before a structure is loaded
  ['colorMode', 'topology', 'scA', 'scB', 'scC', 'scCube', 'boundaryAtoms', 'whiteBg', 'persp', 'savePng',
    'rigidSbu', 'resetDefaults', 'resetZeolite', 'exportOpts', 'importOpts'].forEach((id) => setEnabled(id, loaded));
  document.querySelectorAll('[data-view]').forEach((b) => { b.disabled = !loaded; });
  // the topology view draws neither atoms, H, guests, DOF highlights nor joints
  setEnabled('showDof', loaded && !topo);
  setEnabled('showH', loaded && !topo);
  setEnabled('showGuests', loaded && !topo);
  // joints are only drawn for the undeformed, static structure
  setEnabled('showJoints', loaded && !topo && !moving);
  // MOF-specific options need metals; the zeolite preset needs T-O-T bridges
  setEnabled('rigidSbu', loaded && S.st.has_metal);
  setEnabled('resetDefaults', loaded && S.st.has_carboxylate);
  setEnabled('presetRigid', loaded);
  setEnabled('resetZeolite', loaded && S.st.has_tot);
  // the SBU cutoff matters only while at least one SBU is rigid
  const box = $('rigidSbu');
  setEnabled('sbuCutoff', loaded && S.st.has_metal && (box.checked || box.indeterminate));
  // nothing to follow without a cell-changing infinitesimal mechanism
  setEnabled('followModes', !!S.model && S.model.rigidity.n_macroscopic > 0);
}
updateControlStates();
