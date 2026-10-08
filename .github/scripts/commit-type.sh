#!/usr/bin/env bash
# Decide whether a CI run must be "release strict", i.e. whether R CMD check
# NOTEs (not only ERRORs and WARNINGs) make the run fail.
#
# The run is strict when the pull request title or any commit subject in the
# range being tested is a conventional commit of type `feat` or `fix`, or is a
# breaking change (`type!:` or a `BREAKING CHANGE` footer), and always for
# manual runs (workflow_dispatch), which are used to prepare a release.
#
# Inputs (environment): EVENT_NAME, PR_TITLE, RANGE (a git revision range,
# optional). Outputs (GITHUB_OUTPUT): strict=true|false, reason=<text>.
set -euo pipefail

subjects=""
bodies=""
if [ "${EVENT_NAME:-}" = "pull_request" ] && [ -n "${PR_TITLE:-}" ]; then
  subjects="${PR_TITLE}"$'\n'
fi
if [ -n "${RANGE:-}" ] && git rev-list "${RANGE}" >/dev/null 2>&1; then
  subjects+="$(git log --no-merges --format=%s "${RANGE}")"$'\n'
  bodies="$(git log --no-merges --format=%b "${RANGE}")"
else
  subjects+="$(git log -1 --no-merges --format=%s)"$'\n'
  bodies="$(git log -1 --no-merges --format=%b)"
fi

conv='^[a-z]+(\([^)]+\))?!?: '
strict=false
reason="no feat/fix/breaking commit found"

if [ "${EVENT_NAME:-}" = "workflow_dispatch" ]; then
  strict=true
  reason="manual run (release check)"
fi

while IFS= read -r s; do
  [ -z "$s" ] && continue
  if [[ "$s" =~ ^(feat|fix)(\([^\)]+\))?!?:\  ]] || [[ "$s" =~ ^[a-z]+(\([^\)]+\))?!:\  ]]; then
    strict=true
    reason="conventional commit '${s}'"
    break
  fi
  if ! [[ "$s" =~ $conv ]]; then
    echo "::notice title=Not a conventional commit::${s}"
  fi
done <<< "$subjects"

if [ "$strict" = false ] && grep -q '^BREAKING[ -]CHANGE' <<< "$bodies"; then
  strict=true
  reason="BREAKING CHANGE footer"
fi

echo "strict=${strict}"
echo "reason=${reason}"
if [ -n "${GITHUB_OUTPUT:-}" ]; then
  { echo "strict=${strict}"; echo "reason=${reason}"; } >> "$GITHUB_OUTPUT"
fi
