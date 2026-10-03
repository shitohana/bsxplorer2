#!/bin/sh
set -eu

case "$1" in
    pre-push)
        target_ref=${PRE_COMMIT_REMOTE_BRANCH:-}
        ;;
    pre-merge-commit)
        target_ref=$(git symbolic-ref --quiet HEAD) || exit 1
        ;;
    *)
        echo "Unsupported test hook stage: $1" >&2
        exit 1
        ;;
esac

if [ "$target_ref" != refs/heads/master ]; then
    echo "Skipping Rust tests: target is not master."
    exit 0
fi

exec just rs-test-fast
