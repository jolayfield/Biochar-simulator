#!/usr/bin/env bash
#
# Install the repository's git hooks.
#
# Git hooks live in .git/hooks/, which is not tracked, so every clone has to run
# this once. Uses core.hooksPath so the hooks stay version-controlled and a
# change here reaches everyone without a reinstall.
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

chmod +x tools/hooks/pre-commit tools/hooks/check-secrets.sh
git config core.hooksPath tools/hooks

echo "Installed. git will now run tools/hooks/pre-commit on every commit."
echo
echo "  Checks:  secrets, ruff, fast tests (pytest -m 'not slow'), rq index"
echo "  Bypass:  git commit --no-verify"
echo "  Remove:  git config --unset core.hooksPath"
echo
echo "The full suite -- including the tests marked 'slow' -- runs in CI."
