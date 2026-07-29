#!/usr/bin/env bash
#
# Refuse to commit anything that looks like a credential.
#
# Scans the staged diff only, so it costs nothing and cannot be fooled by a
# clean working tree. This is a tripwire, not a security control -- it catches
# the accidental paste, not a determined one.
#
# Exits non-zero when a match is found.
set -uo pipefail

staged=$(git diff --cached --name-only --diff-filter=ACM)
[ -z "$staged" ] && exit 0

# Added lines only ('^+', excluding the +++ header): an existing match should
# not block every future commit that touches the same file.
added=$(git diff --cached --unified=0 -- $staged | grep -E '^\+' | grep -vE '^\+\+\+' || true)
[ -z "$added" ] && exit 0

hits=0
flag() {
  printf '  possible %s\n' "$1"
  printf '%s\n' "$2" | head -3 | sed 's/^/    /'
  hits=1
}

# Private key blocks.
m=$(printf '%s\n' "$added" | grep -E 'BEGIN (RSA |EC |OPENSSH |PGP )?PRIVATE KEY' || true)
[ -n "$m" ] && flag "private key" "$m"

# Provider-shaped tokens: AWS, GitHub, Slack, Google, generic bearer JWTs.
m=$(printf '%s\n' "$added" | grep -E '(AKIA[0-9A-Z]{16}|gh[pousr]_[A-Za-z0-9]{20,}|xox[abprs]-[A-Za-z0-9-]{10,}|AIza[0-9A-Za-z_-]{35}|eyJ[A-Za-z0-9_-]{10,}\.[A-Za-z0-9_-]{10,}\.)' || true)
[ -n "$m" ] && flag "API token" "$m"

# Assignments to secret-ish names with a non-trivial literal value. Excludes
# obvious placeholders and environment lookups, which are the common false
# positives in config templates and tests.
m=$(printf '%s\n' "$added" \
  | grep -iE '(api[_-]?key|secret|passwd|password|token|credential)[[:space:]]*[:=][[:space:]]*["'"'"'][^"'"'"']{12,}' \
  | grep -viE '(example|placeholder|dummy|changeme|your[_-]?|xxx+|\.\.\.|<[^>]+>|getenv|environ|os\.env)' || true)
[ -n "$m" ] && flag "hardcoded credential" "$m"

exit $hits
