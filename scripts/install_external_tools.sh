#!/usr/bin/env bash
set -uo pipefail

usage() {
  cat <<'EOF'
Install DockMate-VS external tools into the active conda environment.

Usage:
  scripts/install_external_tools.sh [options]

Options:
  --with-pymol          Install open-source PyMOL from conda-forge.
  --package-manager CMD Use this conda-compatible package manager.
  --allow-classic-conda Allow Conda's slower classic solver.
  --no-isolate-engines  Install docking engines into the active environment
                        instead of isolated ones. Vina and Smina have
                        incompatible Boost requirements and cannot coexist,
                        so this will usually fail; kept for troubleshooting.
  --engine-env-dir DIR  Parent directory for isolated engine environments.
                        Defaults to the directory holding the active
                        environment.
  --vina-bin PATH       Register an existing Vina executable instead of installing it.
  --smina-bin PATH      Register an existing Smina executable instead of installing it.
  --rdock-root PATH     Register an existing rDock installation root.
  --pymol-bin PATH      Register an existing PyMOL executable.
  --ligplus-root PATH   Register a separately licensed LigPlot+ installation.
  --dry-run             Print actions without changing the environment.
  -h, --help            Show this help text.

Docking engines (Vina, Smina, rDock) are installed into separate conda
environments by default, because their native dependencies conflict with one
another. Their executables are symlinked into the active environment, so they
remain on PATH exactly as if they had been installed directly.

The script does not edit shell startup files. Conda places the environment's
bin directory on PATH whenever the environment is activated. Manual binaries
are linked into that directory, and activation hooks provide tool-specific
environment variables.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

note() {
  echo "[dockmate-tools] $*"
}

shell_quote() {
  printf '%q' "$1"
}

dry_run=0
with_pymol=0
allow_classic_conda=0
isolate_engines=1
engine_env_dir=""
package_manager_arg=""
vina_bin=""
smina_bin=""
rdock_root=""
pymol_bin=""
ligplus_root=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-pymol)
      with_pymol=1
      shift
      ;;
    --allow-classic-conda)
      allow_classic_conda=1
      shift
      ;;
    --no-isolate-engines)
      isolate_engines=0
      shift
      ;;
    --package-manager|--vina-bin|--smina-bin|--rdock-root|--pymol-bin|--ligplus-root|--engine-env-dir)
      [[ $# -ge 2 ]] || die "$1 requires a path"
      case "$1" in
        --package-manager) package_manager_arg="$2" ;;
        --engine-env-dir) engine_env_dir="$2" ;;
        --vina-bin) vina_bin="$2" ;;
        --smina-bin) smina_bin="$2" ;;
        --rdock-root) rdock_root="$2" ;;
        --pymol-bin) pymol_bin="$2" ;;
        --ligplus-root) ligplus_root="$2" ;;
      esac
      shift 2
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown option: $1"
      ;;
  esac
done

[[ -n "${CONDA_PREFIX:-}" ]] || die "activate the dockmate-vs conda environment first"
[[ "${CONDA_DEFAULT_ENV:-}" != "base" ]] || die "refusing to install external tools into the base environment"
[[ -x "$CONDA_PREFIX/bin/python" ]] || die "Python was not found under $CONDA_PREFIX/bin"

case "$(uname -s)" in
  Darwin|Linux) ;;
  *) die "this installer supports macOS and Linux; use Docker or install tools manually on this platform" ;;
esac

export PATH="$CONDA_PREFIX/bin:$PATH"
note "target environment: $CONDA_PREFIX"

resolve_command() {
  local requested="$1"
  if [[ "$requested" == */* ]]; then
    [[ -x "$requested" ]] || return 1
    printf '%s\n' "$requested"
  else
    command -v "$requested" 2>/dev/null
  fi
}

if [[ -n "$package_manager_arg" ]]; then
  package_manager="$(resolve_command "$package_manager_arg")" || \
    die "package manager is not executable: $package_manager_arg"
elif command -v micromamba >/dev/null 2>&1; then
  package_manager="$(command -v micromamba)"
elif command -v mamba >/dev/null 2>&1; then
  package_manager="$(command -v mamba)"
elif [[ -n "${CONDA_EXE:-}" && -x "$CONDA_EXE" ]]; then
  package_manager="$CONDA_EXE"
else
  package_manager="$(command -v conda || true)"
fi
[[ -x "$package_manager" ]] || die "no conda-compatible package manager was found"
note "package manager: $package_manager"

# Validate supplied locations before starting a potentially expensive solve.
[[ -z "$vina_bin" || -f "$vina_bin" ]] || die "Vina executable was not found: $vina_bin"
[[ -z "$smina_bin" || -f "$smina_bin" ]] || die "Smina executable was not found: $smina_bin"
[[ -z "$pymol_bin" || -f "$pymol_bin" ]] || die "PyMOL executable was not found: $pymol_bin"
[[ -z "$rdock_root" || -d "$rdock_root" ]] || die "rDock root was not found: $rdock_root"
[[ -z "$ligplus_root" || -d "$ligplus_root" ]] || die "LigPlot+ root was not found: $ligplus_root"

command_missing() {
  ! command -v "$1" >/dev/null 2>&1
}

env_command_missing() {
  [[ ! -x "$CONDA_PREFIX/bin/$1" ]]
}

import_missing() {
  ! "$CONDA_PREFIX/bin/python" -c "import $1" >/dev/null 2>&1
}

# Shared packages go into the active environment. Docking engines are queued
# for isolated environments: Vina requires libboost 1.86 while Smina caps at
# libboost 1.82, so a single prefix can never satisfy both.
package_specs=()
engine_queue=()
failed_tools=()

[[ -n "$engine_env_dir" ]] || engine_env_dir="$(dirname "$CONDA_PREFIX")"

queue_engine() {
  local name="$1" spec="$2" binaries="$3"
  if [[ $isolate_engines -eq 1 ]]; then
    engine_queue+=("$name|$spec|$binaries")
  else
    package_specs+=("$spec")
  fi
}

if [[ -z "$vina_bin" ]] && env_command_missing vina; then
  queue_engine vina "vina=1.2.7" "vina"
fi
if [[ -z "$smina_bin" ]] && env_command_missing smina; then
  queue_engine smina "smina=2020.12.10" "smina"
fi
if env_command_missing obabel; then
  package_specs+=("openbabel")
fi
if import_missing openmm; then
  package_specs+=("openmm")
fi
if import_missing pdbfixer; then
  package_specs+=("pdbfixer")
fi
if import_missing propka; then
  package_specs+=("propka")
fi
if env_command_missing reduce; then
  package_specs+=("reduce")
fi
if env_command_missing mk_prepare_receptor.py; then
  package_specs+=("meeko")
fi
if env_command_missing fpocket; then
  package_specs+=("fpocket")
fi
if [[ -z "$rdock_root" ]] && env_command_missing rbdock; then
  if [[ "$(uname -s)" == "Linux" ]]; then
    queue_engine rdock "rdock=24.04.204_legacy" "rbdock,rbcavity"
  else
    note "rDock has no current conda package for macOS; use --rdock-root or Docker."
  fi
fi
if [[ $with_pymol -eq 1 ]] && [[ -z "$pymol_bin" ]] && env_command_missing pymol; then
  package_specs+=("pymol-open-source")
fi
if [[ -n "$ligplus_root" ]]; then
  command_missing java && package_specs+=("openjdk=17")
  command_missing gs && package_specs+=("ghostscript")
fi

if [[ ${#package_specs[@]} -gt 0 || ${#engine_queue[@]} -gt 0 ]]; then
  if [[ "$(basename "$package_manager")" == "conda" ]]; then
    solver_setting="$("$package_manager" config --show solver 2>/dev/null || true)"
    if [[ "$solver_setting" != *"libmamba"* ]]; then
      if [[ $allow_classic_conda -eq 0 ]]; then
        die "Conda is not configured for the libmamba solver. Install Micromamba or Mamba, configure Conda to use libmamba, or rerun with --allow-classic-conda."
      fi
      note "WARNING: using Conda's classic solver; dependency resolution may take a long time."
    fi
  fi
fi

conda_install_into() {
  local prefix="$1"
  shift
  "$package_manager" install --yes --prefix "$prefix" --override-channels \
    --channel conda-forge --channel bioconda "$@"
}

if [[ ${#package_specs[@]} -gt 0 ]]; then
  note "shared packages: ${package_specs[*]}"
  if [[ $dry_run -eq 1 ]]; then
    printf 'DRY RUN:'
    printf ' %q' "$package_manager" install --yes --prefix "$CONDA_PREFIX" \
      --override-channels --channel conda-forge --channel bioconda \
      "rdkit<2026" "${package_specs[@]}"
    printf '\n'
  elif ! conda_install_into "$CONDA_PREFIX" "rdkit<2026" "${package_specs[@]}"; then
    # One unsatisfiable spec must not cost the user every other tool.
    note "WARNING: batch install failed; retrying each package individually"
    for spec in "${package_specs[@]}"; do
      if conda_install_into "$CONDA_PREFIX" "rdkit<2026" "$spec"; then
        note "installed $spec"
      else
        note "WARNING: could not install $spec"
        failed_tools+=("$spec")
      fi
    done
  fi
else
  note "all shared conda-managed tools are already available"
fi

absolute_file() {
  local source_path="$1"
  local source_dir
  [[ -f "$source_path" ]] || return 1
  source_dir="$(cd "$(dirname "$source_path")" && pwd -P)" || return 1
  printf '%s/%s\n' "$source_dir" "$(basename "$source_path")"
}

link_binary() {
  local source_path="$1"
  local target_name="$2"
  local resolved target
  resolved="$(absolute_file "$source_path")" || die "$source_path is not a file"
  [[ -x "$resolved" ]] || die "$resolved is not executable"
  target="$CONDA_PREFIX/bin/$target_name"

  if [[ -e "$target" && ! -L "$target" ]]; then
    die "$target already exists and is not a symlink"
  fi
  if [[ $dry_run -eq 1 ]]; then
    note "DRY RUN: ln -sfn $(shell_quote "$resolved") $(shell_quote "$target")"
  else
    ln -sfn "$resolved" "$target"
  fi
  note "$target_name -> $resolved"
}

[[ -z "$vina_bin" ]] || link_binary "$vina_bin" vina
[[ -z "$smina_bin" ]] || link_binary "$smina_bin" smina
[[ -z "$pymol_bin" ]] || link_binary "$pymol_bin" pymol

engine_env_prefix() {
  printf '%s/dockmate-%s\n' "$engine_env_dir" "$1"
}

# Each engine gets its own prefix so that conflicting native dependencies never
# have to be reconciled, then its executables are linked into the active
# environment's bin directory, which conda already places on PATH.
install_engine() {
  local name="$1" spec="$2" binaries="$3"
  local prefix first_binary binary
  prefix="$(engine_env_prefix "$name")"
  first_binary="${binaries%%,*}"

  if [[ -x "$prefix/bin/$first_binary" ]]; then
    note "reusing isolated environment for $name: $prefix"
  else
    note "installing $name into an isolated environment: $prefix"
    if [[ $dry_run -eq 1 ]]; then
      printf 'DRY RUN:'
      printf ' %q' "$package_manager" create --yes --prefix "$prefix" \
        --override-channels --channel conda-forge --channel bioconda "$spec"
      printf '\n'
      return 0
    fi
    if ! "$package_manager" create --yes --prefix "$prefix" --override-channels \
      --channel conda-forge --channel bioconda "$spec"; then
      note "WARNING: could not install $name; continuing without it"
      failed_tools+=("$name")
      return 1
    fi
  fi

  if [[ "$name" == "rdock" ]]; then
    # configure_rdock links rbdock/rbcavity and exports RBT_ROOT.
    rdock_root="$prefix"
    return 0
  fi

  local IFS=,
  for binary in $binaries; do
    if [[ -x "$prefix/bin/$binary" ]]; then
      link_binary "$prefix/bin/$binary" "$binary"
    else
      note "WARNING: $binary was not found in $prefix"
      failed_tools+=("$binary")
    fi
  done
}

for engine_entry in ${engine_queue[@]+"${engine_queue[@]}"}; do
  IFS='|' read -r engine_name engine_spec engine_binaries <<<"$engine_entry"
  install_engine "$engine_name" "$engine_spec" "$engine_binaries" || true
done

write_activation_pair() {
  local name="$1"
  local variable="$2"
  local value="$3"
  local activate_dir="$CONDA_PREFIX/etc/conda/activate.d"
  local deactivate_dir="$CONDA_PREFIX/etc/conda/deactivate.d"
  local activate_file="$activate_dir/dockmate-$name.sh"
  local deactivate_file="$deactivate_dir/dockmate-$name.sh"
  local quoted_value
  quoted_value="$(shell_quote "$value")"

  if [[ $dry_run -eq 1 ]]; then
    note "DRY RUN: configure $variable=$value on conda activation"
    return
  fi
  mkdir -p "$activate_dir" "$deactivate_dir"
  cat >"$activate_file" <<EOF
export DOCKMATE_PREVIOUS_${variable}="\${${variable}-}"
export ${variable}=${quoted_value}
EOF
  cat >"$deactivate_file" <<EOF
if [[ -n "\${DOCKMATE_PREVIOUS_${variable}:-}" ]]; then
  export ${variable}="\$DOCKMATE_PREVIOUS_${variable}"
else
  unset ${variable}
fi
unset DOCKMATE_PREVIOUS_${variable}
EOF
}

configure_rdock() {
  local root="$1"
  local dock_prm data_dir
  [[ -x "$root/bin/rbdock" ]] || die "rbdock was not found under $root/bin"
  [[ -x "$root/bin/rbcavity" ]] || die "rbcavity was not found under $root/bin"
  if [[ "$root" != "$CONDA_PREFIX" ]]; then
    link_binary "$root/bin/rbdock" rbdock
    link_binary "$root/bin/rbcavity" rbcavity
  fi

  if [[ ! -e "$root/data/scripts/dock.prm" && -d "$root/share" ]]; then
    dock_prm="$(find "$root/share" -path '*/data/scripts/dock.prm' -print -quit 2>/dev/null)"
    if [[ -n "$dock_prm" ]]; then
      data_dir="${dock_prm%/scripts/dock.prm}"
      if [[ $dry_run -eq 1 ]]; then
        note "DRY RUN: link rDock data directory $data_dir to $root/data"
      elif [[ ! -e "$root/data" ]]; then
        ln -s "$data_dir" "$root/data"
      fi
    fi
  fi
  write_activation_pair rdock RBT_ROOT "$root"
  note "rDock root: $root"
}

if [[ -n "$rdock_root" ]]; then
  rdock_root="$(cd "$rdock_root" 2>/dev/null && pwd -P)" || die "invalid rDock root"
  configure_rdock "$rdock_root"
elif [[ -x "$CONDA_PREFIX/bin/rbdock" ]]; then
  configure_rdock "$CONDA_PREFIX"
fi

configure_ligplus() {
  local root="$1"
  local ligplot=""
  local candidate
  root="$(cd "$root" 2>/dev/null && pwd -P)" || die "invalid LigPlot+ root"
  [[ -f "$root/LigPlus.jar" ]] || die "LigPlus.jar was not found under $root"

  if [[ "$(uname -s)" == "Darwin" ]]; then
    for candidate in "$root/lib/exe_mac64/ligplot" "$root/lib/exe_mac/ligplot"; do
      [[ -x "$candidate" ]] && ligplot="$candidate" && break
    done
  else
    for candidate in "$root/lib/exe_linux64/ligplot" "$root/lib/exe_linux/ligplot"; do
      [[ -x "$candidate" ]] && ligplot="$candidate" && break
    done
  fi
  [[ -n "$ligplot" ]] || die "LigPlot executable was not found under $root/lib"

  link_binary "$ligplot" ligplot
  write_activation_pair ligplus-root LIGPLUS_ROOT "$root"
  write_activation_pair ligplus-bin LIGPLOT_BIN "$ligplot"
  write_activation_pair ligplus-jar LIGPLUS_JAR "$root/LigPlus.jar"
  if [[ -f "$root/lib/data/components.cif" ]]; then
    write_activation_pair ligplus-dictionary HET_GROUP_DICTIONARY "$root/lib/data/components.cif"
  fi
  note "LigPlot+ root: $root"
}

if [[ -n "$ligplus_root" ]]; then
  configure_ligplus "$ligplus_root"
fi

if [[ $dry_run -eq 1 ]]; then
  note "dry run complete; no files or environments were changed"
  exit 0
fi

echo
note "executable summary"
for command_name in vina smina obabel reduce mk_prepare_receptor.py fpocket rbdock pymol ligplot; do
  location="$(command -v "$command_name" 2>/dev/null || true)"
  if [[ -n "$location" ]]; then
    printf '  %-24s %s\n' "$command_name" "$location"
  else
    printf '  %-24s %s\n' "$command_name" "not found"
  fi
done

if [[ ${#failed_tools[@]} -gt 0 ]]; then
  echo
  note "the following could not be installed: ${failed_tools[*]}"
  note "supply them with --vina-bin/--smina-bin/--rdock-root, or use the Docker image."
fi

if command_missing vina && command_missing smina; then
  die "neither Vina nor Smina is available"
fi
if command_missing obabel; then
  die "Open Babel (obabel) is unavailable"
fi

note "Conda activation places $CONDA_PREFIX/bin on PATH."
note "Reactivate the environment before using newly written activation hooks."#!/usr/bin/env bash
set -uo pipefail

usage() {
  cat <<'EOF'
Install DockMate-VS external tools into the active conda environment.

Usage:
  scripts/install_external_tools.sh [options]

Options:
  --with-pymol          Install open-source PyMOL from conda-forge.
  --package-manager CMD Use this conda-compatible package manager.
  --allow-classic-conda Allow Conda's slower classic solver.
  --no-isolate-engines  Install docking engines into the active environment
                        instead of isolated ones. Vina and Smina have
                        incompatible Boost requirements and cannot coexist,
                        so this will usually fail; kept for troubleshooting.
  --engine-env-dir DIR  Parent directory for isolated engine environments.
                        Defaults to the directory holding the active
                        environment.
  --vina-bin PATH       Register an existing Vina executable instead of installing it.
  --smina-bin PATH      Register an existing Smina executable instead of installing it.
  --rdock-root PATH     Register an existing rDock installation root.
  --pymol-bin PATH      Register an existing PyMOL executable.
  --ligplus-root PATH   Register a separately licensed LigPlot+ installation.
  --dry-run             Print actions without changing the environment.
  -h, --help            Show this help text.

Docking engines (Vina, Smina, rDock) are installed into separate conda
environments by default, because their native dependencies conflict with one
another. Their executables are symlinked into the active environment, so they
remain on PATH exactly as if they had been installed directly.

The script does not edit shell startup files. Conda places the environment's
bin directory on PATH whenever the environment is activated. Manual binaries
are linked into that directory, and activation hooks provide tool-specific
environment variables.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

note() {
  echo "[dockmate-tools] $*"
}

shell_quote() {
  printf '%q' "$1"
}

dry_run=0
with_pymol=0
allow_classic_conda=0
isolate_engines=1
engine_env_dir=""
package_manager_arg=""
vina_bin=""
smina_bin=""
rdock_root=""
pymol_bin=""
ligplus_root=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-pymol)
      with_pymol=1
      shift
      ;;
    --allow-classic-conda)
      allow_classic_conda=1
      shift
      ;;
    --no-isolate-engines)
      isolate_engines=0
      shift
      ;;
    --package-manager|--vina-bin|--smina-bin|--rdock-root|--pymol-bin|--ligplus-root|--engine-env-dir)
      [[ $# -ge 2 ]] || die "$1 requires a path"
      case "$1" in
        --package-manager) package_manager_arg="$2" ;;
        --engine-env-dir) engine_env_dir="$2" ;;
        --vina-bin) vina_bin="$2" ;;
        --smina-bin) smina_bin="$2" ;;
        --rdock-root) rdock_root="$2" ;;
        --pymol-bin) pymol_bin="$2" ;;
        --ligplus-root) ligplus_root="$2" ;;
      esac
      shift 2
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown option: $1"
      ;;
  esac
done

[[ -n "${CONDA_PREFIX:-}" ]] || die "activate the dockmate-vs conda environment first"
[[ "${CONDA_DEFAULT_ENV:-}" != "base" ]] || die "refusing to install external tools into the base environment"
[[ -x "$CONDA_PREFIX/bin/python" ]] || die "Python was not found under $CONDA_PREFIX/bin"

case "$(uname -s)" in
  Darwin|Linux) ;;
  *) die "this installer supports macOS and Linux; use Docker or install tools manually on this platform" ;;
esac

export PATH="$CONDA_PREFIX/bin:$PATH"
note "target environment: $CONDA_PREFIX"

resolve_command() {
  local requested="$1"
  if [[ "$requested" == */* ]]; then
    [[ -x "$requested" ]] || return 1
    printf '%s\n' "$requested"
  else
    command -v "$requested" 2>/dev/null
  fi
}

if [[ -n "$package_manager_arg" ]]; then
  package_manager="$(resolve_command "$package_manager_arg")" || \
    die "package manager is not executable: $package_manager_arg"
elif command -v micromamba >/dev/null 2>&1; then
  package_manager="$(command -v micromamba)"
elif command -v mamba >/dev/null 2>&1; then
  package_manager="$(command -v mamba)"
elif [[ -n "${CONDA_EXE:-}" && -x "$CONDA_EXE" ]]; then
  package_manager="$CONDA_EXE"
else
  package_manager="$(command -v conda || true)"
fi
[[ -x "$package_manager" ]] || die "no conda-compatible package manager was found"
note "package manager: $package_manager"

# Validate supplied locations before starting a potentially expensive solve.
[[ -z "$vina_bin" || -f "$vina_bin" ]] || die "Vina executable was not found: $vina_bin"
[[ -z "$smina_bin" || -f "$smina_bin" ]] || die "Smina executable was not found: $smina_bin"
[[ -z "$pymol_bin" || -f "$pymol_bin" ]] || die "PyMOL executable was not found: $pymol_bin"
[[ -z "$rdock_root" || -d "$rdock_root" ]] || die "rDock root was not found: $rdock_root"
[[ -z "$ligplus_root" || -d "$ligplus_root" ]] || die "LigPlot+ root was not found: $ligplus_root"

command_missing() {
  ! command -v "$1" >/dev/null 2>&1
}

env_command_missing() {
  [[ ! -x "$CONDA_PREFIX/bin/$1" ]]
}

import_missing() {
  ! "$CONDA_PREFIX/bin/python" -c "import $1" >/dev/null 2>&1
}

# Shared packages go into the active environment. Docking engines are queued
# for isolated environments: Vina requires libboost 1.86 while Smina caps at
# libboost 1.82, so a single prefix can never satisfy both.
package_specs=()
engine_queue=()
failed_tools=()

[[ -n "$engine_env_dir" ]] || engine_env_dir="$(dirname "$CONDA_PREFIX")"

queue_engine() {
  local name="$1" spec="$2" binaries="$3"
  if [[ $isolate_engines -eq 1 ]]; then
    engine_queue+=("$name|$spec|$binaries")
  else
    package_specs+=("$spec")
  fi
}

if [[ -z "$vina_bin" ]] && env_command_missing vina; then
  queue_engine vina "vina=1.2.7" "vina"
fi
if [[ -z "$smina_bin" ]] && env_command_missing smina; then
  queue_engine smina "smina=2020.12.10" "smina"
fi
if env_command_missing obabel; then
  package_specs+=("openbabel")
fi
if import_missing openmm; then
  package_specs+=("openmm")
fi
if import_missing pdbfixer; then
  package_specs+=("pdbfixer")
fi
if import_missing propka; then
  package_specs+=("propka")
fi
if env_command_missing reduce; then
  package_specs+=("reduce")
fi
if env_command_missing mk_prepare_receptor.py; then
  package_specs+=("meeko")
fi
if env_command_missing fpocket; then
  package_specs+=("fpocket")
fi
if [[ -z "$rdock_root" ]] && env_command_missing rbdock; then
  if [[ "$(uname -s)" == "Linux" ]]; then
    queue_engine rdock "rdock=24.04.204_legacy" "rbdock,rbcavity"
  else
    note "rDock has no current conda package for macOS; use --rdock-root or Docker."
  fi
fi
if [[ $with_pymol -eq 1 ]] && [[ -z "$pymol_bin" ]] && env_command_missing pymol; then
  package_specs+=("pymol-open-source")
fi
if [[ -n "$ligplus_root" ]]; then
  command_missing java && package_specs+=("openjdk=17")
  command_missing gs && package_specs+=("ghostscript")
fi

if [[ ${#package_specs[@]} -gt 0 || ${#engine_queue[@]} -gt 0 ]]; then
  if [[ "$(basename "$package_manager")" == "conda" ]]; then
    solver_setting="$("$package_manager" config --show solver 2>/dev/null || true)"
    if [[ "$solver_setting" != *"libmamba"* ]]; then
      if [[ $allow_classic_conda -eq 0 ]]; then
        die "Conda is not configured for the libmamba solver. Install Micromamba or Mamba, configure Conda to use libmamba, or rerun with --allow-classic-conda."
      fi
      note "WARNING: using Conda's classic solver; dependency resolution may take a long time."
    fi
  fi
fi

conda_install_into() {
  local prefix="$1"
  shift
  "$package_manager" install --yes --prefix "$prefix" --override-channels \
    --channel conda-forge --channel bioconda "$@"
}

if [[ ${#package_specs[@]} -gt 0 ]]; then
  note "shared packages: ${package_specs[*]}"
  if [[ $dry_run -eq 1 ]]; then
    printf 'DRY RUN:'
    printf ' %q' "$package_manager" install --yes --prefix "$CONDA_PREFIX" \
      --override-channels --channel conda-forge --channel bioconda \
      "rdkit<2026" "${package_specs[@]}"
    printf '\n'
  elif ! conda_install_into "$CONDA_PREFIX" "rdkit<2026" "${package_specs[@]}"; then
    # One unsatisfiable spec must not cost the user every other tool.
    note "WARNING: batch install failed; retrying each package individually"
    for spec in "${package_specs[@]}"; do
      if conda_install_into "$CONDA_PREFIX" "rdkit<2026" "$spec"; then
        note "installed $spec"
      else
        note "WARNING: could not install $spec"
        failed_tools+=("$spec")
      fi
    done
  fi
else
  note "all shared conda-managed tools are already available"
fi

absolute_file() {
  local source_path="$1"
  local source_dir
  [[ -f "$source_path" ]] || return 1
  source_dir="$(cd "$(dirname "$source_path")" && pwd -P)" || return 1
  printf '%s/%s\n' "$source_dir" "$(basename "$source_path")"
}

link_binary() {
  local source_path="$1"
  local target_name="$2"
  local resolved target
  resolved="$(absolute_file "$source_path")" || die "$source_path is not a file"
  [[ -x "$resolved" ]] || die "$resolved is not executable"
  target="$CONDA_PREFIX/bin/$target_name"

  if [[ -e "$target" && ! -L "$target" ]]; then
    die "$target already exists and is not a symlink"
  fi
  if [[ $dry_run -eq 1 ]]; then
    note "DRY RUN: ln -sfn $(shell_quote "$resolved") $(shell_quote "$target")"
  else
    ln -sfn "$resolved" "$target"
  fi
  note "$target_name -> $resolved"
}

[[ -z "$vina_bin" ]] || link_binary "$vina_bin" vina
[[ -z "$smina_bin" ]] || link_binary "$smina_bin" smina
[[ -z "$pymol_bin" ]] || link_binary "$pymol_bin" pymol

engine_env_prefix() {
  printf '%s/dockmate-%s\n' "$engine_env_dir" "$1"
}

# Each engine gets its own prefix so that conflicting native dependencies never
# have to be reconciled, then its executables are linked into the active
# environment's bin directory, which conda already places on PATH.
install_engine() {
  local name="$1" spec="$2" binaries="$3"
  local prefix first_binary binary
  prefix="$(engine_env_prefix "$name")"
  first_binary="${binaries%%,*}"

  if [[ -x "$prefix/bin/$first_binary" ]]; then
    note "reusing isolated environment for $name: $prefix"
  else
    note "installing $name into an isolated environment: $prefix"
    if [[ $dry_run -eq 1 ]]; then
      printf 'DRY RUN:'
      printf ' %q' "$package_manager" create --yes --prefix "$prefix" \
        --override-channels --channel conda-forge --channel bioconda "$spec"
      printf '\n'
      return 0
    fi
    if ! "$package_manager" create --yes --prefix "$prefix" --override-channels \
      --channel conda-forge --channel bioconda "$spec"; then
      note "WARNING: could not install $name; continuing without it"
      failed_tools+=("$name")
      return 1
    fi
  fi

  if [[ "$name" == "rdock" ]]; then
    # configure_rdock links rbdock/rbcavity and exports RBT_ROOT.
    rdock_root="$prefix"
    return 0
  fi

  local IFS=,
  for binary in $binaries; do
    if [[ -x "$prefix/bin/$binary" ]]; then
      link_binary "$prefix/bin/$binary" "$binary"
    else
      note "WARNING: $binary was not found in $prefix"
      failed_tools+=("$binary")
    fi
  done
}

for engine_entry in ${engine_queue[@]+"${engine_queue[@]}"}; do
  IFS='|' read -r engine_name engine_spec engine_binaries <<<"$engine_entry"
  install_engine "$engine_name" "$engine_spec" "$engine_binaries" || true
done

write_activation_pair() {
  local name="$1"
  local variable="$2"
  local value="$3"
  local activate_dir="$CONDA_PREFIX/etc/conda/activate.d"
  local deactivate_dir="$CONDA_PREFIX/etc/conda/deactivate.d"
  local activate_file="$activate_dir/dockmate-$name.sh"
  local deactivate_file="$deactivate_dir/dockmate-$name.sh"
  local quoted_value
  quoted_value="$(shell_quote "$value")"

  if [[ $dry_run -eq 1 ]]; then
    note "DRY RUN: configure $variable=$value on conda activation"
    return
  fi
  mkdir -p "$activate_dir" "$deactivate_dir"
  cat >"$activate_file" <<EOF
export DOCKMATE_PREVIOUS_${variable}="\${${variable}-}"
export ${variable}=${quoted_value}
EOF
  cat >"$deactivate_file" <<EOF
if [[ -n "\${DOCKMATE_PREVIOUS_${variable}:-}" ]]; then
  export ${variable}="\$DOCKMATE_PREVIOUS_${variable}"
else
  unset ${variable}
fi
unset DOCKMATE_PREVIOUS_${variable}
EOF
}

configure_rdock() {
  local root="$1"
  local dock_prm data_dir
  [[ -x "$root/bin/rbdock" ]] || die "rbdock was not found under $root/bin"
  [[ -x "$root/bin/rbcavity" ]] || die "rbcavity was not found under $root/bin"
  if [[ "$root" != "$CONDA_PREFIX" ]]; then
    link_binary "$root/bin/rbdock" rbdock
    link_binary "$root/bin/rbcavity" rbcavity
  fi

  if [[ ! -e "$root/data/scripts/dock.prm" && -d "$root/share" ]]; then
    dock_prm="$(find "$root/share" -path '*/data/scripts/dock.prm' -print -quit 2>/dev/null)"
    if [[ -n "$dock_prm" ]]; then
      data_dir="${dock_prm%/scripts/dock.prm}"
      if [[ $dry_run -eq 1 ]]; then
        note "DRY RUN: link rDock data directory $data_dir to $root/data"
      elif [[ ! -e "$root/data" ]]; then
        ln -s "$data_dir" "$root/data"
      fi
    fi
  fi
  write_activation_pair rdock RBT_ROOT "$root"
  note "rDock root: $root"
}

if [[ -n "$rdock_root" ]]; then
  rdock_root="$(cd "$rdock_root" 2>/dev/null && pwd -P)" || die "invalid rDock root"
  configure_rdock "$rdock_root"
elif [[ -x "$CONDA_PREFIX/bin/rbdock" ]]; then
  configure_rdock "$CONDA_PREFIX"
fi

configure_ligplus() {
  local root="$1"
  local ligplot=""
  local candidate
  root="$(cd "$root" 2>/dev/null && pwd -P)" || die "invalid LigPlot+ root"
  [[ -f "$root/LigPlus.jar" ]] || die "LigPlus.jar was not found under $root"

  if [[ "$(uname -s)" == "Darwin" ]]; then
    for candidate in "$root/lib/exe_mac64/ligplot" "$root/lib/exe_mac/ligplot"; do
      [[ -x "$candidate" ]] && ligplot="$candidate" && break
    done
  else
    for candidate in "$root/lib/exe_linux64/ligplot" "$root/lib/exe_linux/ligplot"; do
      [[ -x "$candidate" ]] && ligplot="$candidate" && break
    done
  fi
  [[ -n "$ligplot" ]] || die "LigPlot executable was not found under $root/lib"

  link_binary "$ligplot" ligplot
  write_activation_pair ligplus-root LIGPLUS_ROOT "$root"
  write_activation_pair ligplus-bin LIGPLOT_BIN "$ligplot"
  write_activation_pair ligplus-jar LIGPLUS_JAR "$root/LigPlus.jar"
  if [[ -f "$root/lib/data/components.cif" ]]; then
    write_activation_pair ligplus-dictionary HET_GROUP_DICTIONARY "$root/lib/data/components.cif"
  fi
  note "LigPlot+ root: $root"
}

if [[ -n "$ligplus_root" ]]; then
  configure_ligplus "$ligplus_root"
fi

if [[ $dry_run -eq 1 ]]; then
  note "dry run complete; no files or environments were changed"
  exit 0
fi

echo
note "executable summary"
for command_name in vina smina obabel reduce mk_prepare_receptor.py fpocket rbdock pymol ligplot; do
  location="$(command -v "$command_name" 2>/dev/null || true)"
  if [[ -n "$location" ]]; then
    printf '  %-24s %s\n' "$command_name" "$location"
  else
    printf '  %-24s %s\n' "$command_name" "not found"
  fi
done

if [[ ${#failed_tools[@]} -gt 0 ]]; then
  echo
  note "the following could not be installed: ${failed_tools[*]}"
  note "supply them with --vina-bin/--smina-bin/--rdock-root, or use the Docker image."
fi

if command_missing vina && command_missing smina; then
  die "neither Vina nor Smina is available"
fi
if command_missing obabel; then
  die "Open Babel (obabel) is unavailable"
fi

note "Conda activation places $CONDA_PREFIX/bin on PATH."
note "Reactivate the environment before using newly written activation hooks."
