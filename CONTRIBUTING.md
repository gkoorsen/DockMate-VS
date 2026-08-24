# Contributing

Contributions that improve correctness, reproducibility, documentation, or
platform support are welcome.

## Report a problem

Open a GitHub issue with:

- the application and docking-engine versions;
- operating system and Python version;
- the smallest non-confidential workbook that reproduces the problem;
- the run manifest, relevant result rows, and complete error message; and
- whether the run was new, resumed, or reopened for reporting.

Do not upload proprietary compounds, credentials, or restricted datasets.

## Develop locally

```bash
git clone https://github.com/gkoorsen/docking_platform_gui.git
cd docking_platform_gui
conda env create -f environment.yml
conda activate docking-platform-gui
python -m pip install -e ".[test]"
python -m pytest docking_platform_gui/tests -q
```

Use a focused branch and keep each pull request limited to one concern. Add a
regression test for behavioural changes. Preserve machine-readable output
fields when possible; document any schema change in `CHANGELOG.md`.

## Pull-request checklist

- Tests pass locally and new behaviour has coverage.
- User-visible changes are documented in `README.md` or `USER_GUIDE.md`.
- No generated campaign outputs, large datasets, or external binaries are
  committed.
- New dependencies have a clear need and compatible license.
- Scientific claims distinguish platform behaviour from docking-engine
  behaviour.

By participating, you agree to follow `CODE_OF_CONDUCT.md`.
