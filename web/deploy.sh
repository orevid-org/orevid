#!/bin/bash
# Publishes the static website (web/public, including the generated data) to the gh-pages branch of the
# repository, as a single commit (old versions are not kept, so the repository does not grow).
#   bash web/deploy.sh                    # remote "origin"
#   DOMAIN=orevid.org bash web/deploy.sh  # also writes the CNAME file for a custom domain
set -e
cd "$(dirname "$0")/.."
REMOTE="${REMOTE:-$(git remote get-url origin)}"
[ -f web/public/data/index.json ] || { echo "web/public/data is missing: run python -m web.build first"; exit 1; }
TMP="$(mktemp -d)"
cp -r web/public/. "$TMP"/
touch "$TMP/.nojekyll"                     # serve files as they are (no Jekyll processing)
[ -n "$DOMAIN" ] && echo "$DOMAIN" > "$TMP/CNAME"
cd "$TMP"
git init -q -b gh-pages
git add -A
git -c user.name="OREVID deploy" -c user.email="deploy@orevid.org" commit -q -m "Website build $(date -u +%Y-%m-%d)"
git push -f "$REMOTE" gh-pages
echo "published $(git ls-files | wc -l) files to gh-pages of $REMOTE"
cd - > /dev/null
rm -rf "$TMP"
