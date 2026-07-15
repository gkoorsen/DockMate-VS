# Setting up the Docking GUI on Windows with WSL

This guide is for **students on a Windows PC**. It walks you from a fresh
machine to a running Docking GUI, using **WSL** (Windows Subsystem for Linux)
so you get a real Linux environment where all the docking tools work.

Follow the steps in order. Commands that run in **Windows PowerShell** and
commands that run **inside Ubuntu (WSL)** are labelled every time — don't mix
them up.

**Time needed:** ~30–60 minutes, most of it waiting for downloads.

---

## What you'll end up with

- Ubuntu running inside Windows (via WSL 2)
- A conda environment called `docking_platform_gui` containing Python, RDKit,
  the docking engines, and every other tool the GUI needs
- The GUI window opening on your Windows desktop

---

## Prerequisites

- **Windows 11**, or **Windows 10 version 22H2** or newer. The graphical window
  relies on **WSLg**, which is built into these versions.
  (Check your version: press `Win + R`, type `winver`, press Enter.)
- About **15 GB** of free disk space.
- **Access to the repository.** If the repository is private, ask the project
  owner to add your GitHub account as a collaborator first (see
  [Step 4](#step-4-get-the-code-inside-ubuntu)).

---

## Step 1 — Install WSL and Ubuntu

**In Windows PowerShell (run as Administrator** — right-click the Start button →
"Terminal (Admin)" or "Windows PowerShell (Admin)"**):**

```powershell
wsl --install
```

This installs WSL 2 and Ubuntu. **Restart your PC** when it asks.

After the restart, an **Ubuntu** terminal window opens automatically (if not,
open the "Ubuntu" app from the Start menu). The first time it runs you'll be
asked to create a **username** and **password** — pick anything you'll
remember. This password is used for `sudo` (administrator) commands later; it's
normal that nothing appears on screen while you type it.

Then make sure WSL is fully up to date. **Back in PowerShell:**

```powershell
wsl --update
wsl --shutdown
```

Re-open the **Ubuntu** terminal. **Everything from here on runs inside Ubuntu**,
unless it says otherwise.

Update the Ubuntu package lists and install the couple of system tools we need:

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y git wget
```

---

## Step 2 — Check the graphical display works

The GUI needs to be able to open a window. On Windows 11 / Windows 10 22H2 this
works out of the box through WSLg. Test it:

```bash
sudo apt install -y x11-apps
xeyes
```

A little pair of googly eyes should pop up on your desktop. If it does,
graphics work — close it (click the window's X, or press `Ctrl + C` in the
terminal) and continue.

> **If no window appears:** run `wsl --update` in PowerShell again, then
> `wsl --shutdown`, reopen Ubuntu and retry. Make sure Windows itself is fully
> updated (Settings → Windows Update). WSLg requires Windows 11 or Windows 10
> 22H2+. See [Troubleshooting](#troubleshooting).

---

## Step 3 — Install Miniconda

Conda installs Python **and** the compiled docking tools (smina, Open Babel,
fpocket, OpenMM, …) that plain `pip` cannot.

**Inside Ubuntu:**

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda3
~/miniconda3/bin/conda init bash
```

Close the Ubuntu terminal and open it again (this activates conda). You should
now see `(base)` at the start of your prompt. Verify:

```bash
conda --version
```

---

## Step 4 — Get the code (inside Ubuntu)

> **Important — where to put it.** Clone into your Linux home folder (`~`), as
> shown below. Do **not** put it under `/mnt/c/...` (your Windows C: drive) —
> it will be much slower and can cause permission errors.

```bash
mkdir -p ~/projects
cd ~/projects
git clone https://github.com/gkoorsen/docking_platform_gui.git
cd docking_platform_gui
```

**Then switch to the branch that has the Linux setup fixes:**

```bash
git checkout claude/linux-install-fixes
```

> **Why a branch?** The install/setup fixes (the conda `environment.yml`,
> packaging, etc.) currently live on the `claude/linux-install-fixes` branch.
> Once that branch is merged into `master`, you can skip this `git checkout`
> and just use `master`.

### If the repository is private

The `git clone` above will ask you to sign in. The simplest way:

1. Install the GitHub CLI and authenticate once — it stores your login for
   all future `git` commands:

   ```bash
   sudo apt install -y gh
   gh auth login
   ```

   Choose **GitHub.com → HTTPS → "Login with a web browser"** and follow the
   prompts. Then run the `git clone` command again.

2. Alternatively, when `git` asks for a password, paste a
   **Personal Access Token** (GitHub → Settings → Developer settings →
   Personal access tokens), not your account password.

---

## Step 5 — Create the environment

Still inside `~/projects/docking_platform_gui`:

```bash
conda env create -f environment.yml
```

This downloads and installs everything, including the package itself. **It can
take 10–20 minutes** — it's normal for the "Solving environment" step to sit
for a while. Leave it running.

When it finishes, activate the environment:

```bash
conda activate docking_platform_gui
```

Your prompt should now start with `(docking_platform_gui)`. **You must run this
`conda activate` command every time you open a new Ubuntu terminal to use the
GUI.**

---

## Step 6 — Launch the GUI

```bash
redock-gui
```

The Docking GUI window should open on your Windows desktop. 🎉

(Equivalently: `python scripts/launch_redock_analysis_gui.py`.)

---

## Using it afterwards

Every time you want to run the GUI again:

1. Open the **Ubuntu** terminal.
2. Run these two lines:

   ```bash
   cd ~/projects/docking_platform_gui
   conda activate docking_platform_gui
   redock-gui
   ```

### Getting the latest code later

```bash
cd ~/projects/docking_platform_gui
git pull
```

If `environment.yml` changed after an update, refresh the environment:

```bash
conda env update -f environment.yml --prune
```

---

## Troubleshooting

**The GUI window never appears / `xeyes` didn't work.**
Your Windows may be too old or WSLg isn't active. In PowerShell run
`wsl --update` then `wsl --shutdown`, reopen Ubuntu, and try again. Confirm you
are on Windows 11 or Windows 10 22H2+ (`winver`). Make sure `echo $DISPLAY`
inside Ubuntu prints something like `:0` — if it's empty, WSLg isn't running.

**`conda: command not found`.**
You skipped reopening the terminal after `conda init`. Close and reopen Ubuntu,
or run `source ~/.bashrc`.

**`conda env create` is extremely slow or fails to "solve".**
Recent Miniconda uses the fast `libmamba` solver by default. If yours is old,
update conda first: `conda update -n base -c defaults conda`. If one specific
tool refuses to install for your system, you can open `environment.yml`, put a
`#` in front of that one line to skip it, and re-run the command — the GUI will
still start (you just won't have that one optional tool).

**`redock-gui: command not found`.**
The environment isn't active. Run `conda activate docking_platform_gui` first.
The `(docking_platform_gui)` prefix must be showing in your prompt.

**Docking says a tool like `smina` or `obabel` is missing.**
Make sure the environment is active (see above). These tools are installed by
`environment.yml` inside the `docking_platform_gui` environment only.

**"Protein preparation requires PDBFixer and OpenMM…" error.**
Those come with `environment.yml`. If you see this, your environment isn't
active, or you built a pip-only setup — activate the conda environment.

**Everything feels slow.**
Make sure you cloned into `~/projects` (Linux side), not `/mnt/c/...` (Windows
side). Move it with `mv /mnt/c/…/docking_platform_gui ~/projects/` if needed.

---

## Quick reference (all the commands in order)

```powershell
# --- Windows PowerShell (Admin) ---
wsl --install          # then restart Windows, create Ubuntu username/password
wsl --update
wsl --shutdown
```

```bash
# --- Ubuntu (WSL) ---
sudo apt update && sudo apt upgrade -y
sudo apt install -y git wget x11-apps
xeyes                                   # test graphics, then close it

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda3
~/miniconda3/bin/conda init bash
# close and reopen the Ubuntu terminal here

mkdir -p ~/projects && cd ~/projects
git clone https://github.com/gkoorsen/docking_platform_gui.git
cd docking_platform_gui
git checkout claude/linux-install-fixes

conda env create -f environment.yml
conda activate docking_platform_gui
redock-gui
```
