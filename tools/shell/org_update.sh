#!/usr/bin/env bash
set -euo pipefail

# Tangle rna-seq.org and export README
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ORG_FILE="${REPO_DIR}/rna-seq.org"

emacsclient --socket-name ~/.emacs.d/server/server \
  --eval "(progn
    (find-file \"${ORG_FILE}\")
    (revert-buffer t t)
    (org-babel-tangle)
    (org-md-export-to-markdown))"
