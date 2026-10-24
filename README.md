# NetsUnhinged

A vibe-coded web app that predicts whether a framework material (MOFs, ZIFs, zeolites and other nets)
is **potentially flexible** by treating it as a periodic network of rigid bodies joined by hinges
and pivots. It implements and generalises

> L. Sarkisov, R. L. Martin, M. Haranczyk, B. Smit, *On the Flexibility of Metal–Organic
> Frameworks*, J. Am. Chem. Soc. **2014**, 136, 2228, [doi:10.1021/ja411673b](https://doi.org/10.1021/ja411673b).

What it does:

* reads a CIF (symmetry expansion, removal of guests, H and other dangling atoms),
* shows the structure in 3D, coloured by element, rigid body, or degrees of freedom,
* lets you **click atoms and bonds to choose what is flexible**,
* analyses the resulting mechanism twice: once infinitesimally (exact linear algebra) and once with the
  paper's finite deformation test (13 cell-perturbation regimes).

## Installation

NetsUnhinged needs Python ≥ 3.10 and installs `numpy` and `scipy`. three.js (0.160.0, MIT licence)
is bundled, so the viewer works offline and makes no third-party requests.

NetsUnhinged is on [PyPI](https://pypi.org/project/netsunhinged/).

**Without Python, with uv (simplest).** [uv](https://docs.astral.sh/uv/) downloads Python by
itself. Install uv once (macOS / Linux):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

or on Windows (PowerShell):

```bash
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Then start the app at any time with:

```bash
uvx netsunhinged
```

or install the `netsunhinged` command permanently with `uv tool install netsunhinged`. If the
command is not found afterwards, run `uv tool update-shell` and open a new terminal.

**With pipx.** [pipx](https://pipx.pypa.io) installs the app in its own environment and puts the
`netsunhinged` command on your PATH:

```bash
pipx install netsunhinged
```

**With pip**, ideally inside a virtual environment or conda environment:

```bash
pip install netsunhinged
```

**With Docker** (then open http://localhost:8765):

```bash
docker run --rm -p 8765:8765 ghcr.io/jbweiss-chem/netsunhinged
```

**The latest development version** from GitHub:

```bash
pipx install git+https://github.com/JBWeiss-chem/netsunhinged
```

To update: `uv tool upgrade netsunhinged`, `pipx upgrade netsunhinged`, `pip install --upgrade
netsunhinged` or `docker pull` the image. `uvx` reuses the version it downloaded first; run
`uvx netsunhinged@latest` to get the newest release.

## Usage

After installing, start the app with:

```bash
netsunhinged
```

Your browser opens the app at http://127.0.0.1:8765. It runs only on your computer: nothing is
uploaded anywhere. Stop it with Ctrl+C in the terminal. `netsunhinged serve` does the same and takes
options: `--port 9000` uses another port, `--no-browser` does not open a browser window.

If the command is not found (e.g. because the environment you installed it into is not active), use
`python -m netsunhinged` instead.

Analyse CIF files from the command line:

```bash
netsunhinged analyze my.cif other.cif --preset paper
```

`netsunhinged analyze --help` lists all options. `python -m netsunhinged` works as well.

The **Examples** menu lists the literature structures in `netsunhinged/examples/`: open-licence (CC0)
copies from the Crystallography Open Database, named `{trivial name}_COD-{id}.cif` (see
`netsunhinged/examples/README.md`). In a source checkout, files in `examples_local/` (e.g. CSD
downloads, which may not be redistributed) appear locally as "Local only" and are never published.

## Limitations

* This is a screening tool. As the paper notes, it finds *potential* modes. Energetics, sterics and
  guest effects need force-field or DFT follow-up.
* The model is all-or-nothing: an angle or torsion is either fixed or completely free. A mechanism
  or a "flexible" regime therefore means that the deformation is possible without distorting any
  rigid body, not that it costs no energy. In the material, the released angles and torsions,
  non-bonded interactions and guests still resist it. The reported E and F measure the geometric
  mismatch at the joints (Å²), not an energy.
* Bodies are geometric: atoms have no size, so very large deformations can pass through each other.
  Inspect the deformed frames.
* Bonds come from covalent radii (+0.45 Å). Only the first disorder component is kept, and partially
  occupied sites are kept as full atoms.
* Frameworks with fewer than 3 periodic dimensions (layers, chains) have no bonds, and so no
  constraints, in the remaining directions. The model reports them as flexible there, although
  interlayer or interchain interactions may hold them in place. A warning is shown.

## License

MIT, see `LICENSE`. The bundled three.js is MIT-licensed as well. The example structures come from
the Crystallography Open Database and are in the public domain (CC0); please cite the original
publications listed in `netsunhinged/examples/README.md`.
