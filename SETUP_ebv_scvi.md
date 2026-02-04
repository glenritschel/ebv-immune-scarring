# EBV Immune Scarring – Environment Bootstrap (Clean)

This document recreates the full computational environment for the
`ebv-immune-scarring` project on a fresh WSL Ubuntu install using micromamba.

Environment name: `ebv-scvi`

---

## System prerequisites (one time)

```bash
sudo apt update
sudo apt install -y \
    curl \
    ca-certificates \
    tar \
    bzip2 \
    build-essential \
    g++ \
    make
````

---

## Install micromamba (if not already installed)

```bash
cd ~
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
sudo mv bin/micromamba /usr/local/bin/
micromamba shell init -s bash
source ~/.bashrc
```

Verify:

```bash
micromamba --version
```

---

## Fix broken conda config (if present)

```bash
mv ~/.condarc ~/.condarc.bak
```

---

## Clone repository

```bash
mkdir -p ~/repos
cd ~/repos

git clone git@github.com:glenritschel/ebv-immune-scarring.git
cd ebv-immune-scarring
```

---

## Create EBV environment

```bash
micromamba env create -n ebv-scvi -f environment.yml
micromamba activate ebv-scvi
```

---

## Install pip-only dependencies

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements-pip.txt
```

---

## Register Jupyter kernel (recommended)

```bash
python -m pip install ipykernel
python -m ipykernel install --user --name ebv-scvi --display-name "ebv-scvi"
```

---

## Validate environment

```bash
python -c "import numpy, pandas, scanpy, anndata; print('environment OK')"
```

If scvi-tools is used:

```bash
python -c "import scvi; print('scvi OK')"
```

---

## Reproduction instructions

See:

```bash
REPRODUCE.md
```

for exact figure and analysis commands.

---

## Operational rules (important)

Always shut down WSL before reboot or sleep:

```bash
wsl --shutdown
```

Never unplug the external SSD while Windows is running.
Do not back up live `X:\wsl` VHDX files.
Keep active repos inside WSL (`~/repos`).

---

## Result

* Clean WSL environment on external SSD
* Dedicated `ebv-scvi` micromamba environment
* Reproducible EBV immune scarring pipeline
* Jupyter-ready kernel

```

