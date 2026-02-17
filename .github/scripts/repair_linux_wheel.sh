#!/usr/bin/env bash
set -euo pipefail

wheel_path="${1:?missing wheel path}"
dest_dir="${2:?missing destination dir}"
search_roots="${SEARCH_ROOTS:-/project /opt/rdkit-deps /tmp}"

# 1) locate RDKit libs in likely roots (repo + deps prefix + tmp)
RDKIT_LIB_DIRS="$({
  find ${search_roots} -type f -name 'libRDKit*.so*' -printf '%h\n' 2>/dev/null || true
} | sort -u | tr '\n' ':')"

echo "RDKIT_LIB_DIRS=${RDKIT_LIB_DIRS}"

# 2) make sure SONAME symlinks exist next to versioned files
IFS=':' read -r -a dirs <<< "${RDKIT_LIB_DIRS}"
for d in "${dirs[@]}"; do
  [ -d "$d" ] || continue
  while IFS= read -r f; do
    b="$(basename "$f")"
    soname="$(printf '%s' "$b" | sed -E 's/^(libRDKit[^ ]*\.so\.[0-9]+).*/\1/')"
    if [ -n "$soname" ] && [ ! -e "$d/$soname" ]; then
      ln -s "$b" "$d/$soname" || true
    fi
  done < <(find "$d" -maxdepth 1 -type f -name 'libRDKit*.so.*')
done

# 3) sanity check for common missing lib
if ! (
  echo "$RDKIT_LIB_DIRS" | tr ':' '\n' | while IFS= read -r d; do
    [ -e "$d/libRDKitRDGeometryLib.so.1" ] && echo "FOUND: $d" && exit 0
  done
); then
  echo "ERROR: libRDKitRDGeometryLib.so.1 still not found. Dumping candidates:"
  find ${search_roots} -name 'libRDKitRDGeometryLib.so*' -print || true
  exit 1
fi

export LD_LIBRARY_PATH="${RDKIT_LIB_DIRS}/opt/rdkit-deps/lib:${LD_LIBRARY_PATH:-}"

auditwheel show "$wheel_path"
auditwheel repair -w "$dest_dir" "$wheel_path"
