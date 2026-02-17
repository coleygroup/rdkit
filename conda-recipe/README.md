# Build and install this fork with mamba

This recipe packages the local source tree as a conda package named `rdkit-typer`.

## 1) Create a build environment

```bash
mamba create -n rdkit_pkg -c conda-forge \
  mamba conda-build boa \
  python=3.11
mamba activate rdkit_pkg
```

## 2) Build the package

From the repository root:

```bash
cd /home/bmahjour/test_rdpp/rdkit
conda mambabuild conda-recipe
```

The package artifact is produced in your conda-bld directory (typically `~/miniforge3/conda-bld/linux-64/`).

## 3) Install it with mamba

Option A: install directly from local `conda-bld`:

```bash
mamba install -n myenv -c local rdkit-typer
```

Option B: index a custom folder and install from it:

```bash
mkdir -p /tmp/rdkit-channel/linux-64
cp ~/miniforge3/conda-bld/linux-64/rdkit-typer-*.conda /tmp/rdkit-channel/linux-64/
conda index /tmp/rdkit-channel
mamba install -n myenv -c file:///tmp/rdkit-channel rdkit-typer
```

## 4) Quick check

```bash
python -c "from rdkit import Chem; print(Chem.MolFromSmiles('CCO') is not None)"
```
