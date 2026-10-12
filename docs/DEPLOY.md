# Publishing the website

The website is static: `web/public/` (HTML, JavaScript, CSS, locale files) plus the generated data in
`web/public/data/` (~100 MB, one JSON fragment per 500 diseases). There is no server-side code.

## Build

```bash
python -m web.build        # all diseases with genes or phenotypes (~1 h); needs the built graph (build_all.sh)
python web/serve.py        # check it locally at http://localhost:8796
```

## Publish to GitHub Pages

```bash
bash web/deploy.sh                     # pushes web/public as a single commit to the gh-pages branch
DOMAIN=orevid.org bash web/deploy.sh   # same, with the custom domain (writes CNAME)
```

Then, once, in the repository settings: *Settings → Pages → Source: Deploy from a branch → `gh-pages`,
folder `/ (root)`*. With a custom domain, set it under *Custom domain* and tick *Enforce HTTPS*.

## Custom domain (DNS at the registrar)

| Type | Host | Value |
|---|---|---|
| A | @ | 185.199.108.153 |
| A | @ | 185.199.109.153 |
| A | @ | 185.199.110.153 |
| A | @ | 185.199.111.153 |
| AAAA | @ | 2606:50c0:8000::153, 2606:50c0:8001::153, 2606:50c0:8002::153, 2606:50c0:8003::153 |
| CNAME | www | `<organisation>.github.io` |

Verify the domain for the organisation (*Organisation settings → Pages → Add a domain*, a TXT record) so
that no other repository can claim it.

## Updating

When a source releases new data (e.g. a new Open Targets version), rebuild with `bash build_all.sh`,
check the panel (`docs/PANEL.md`) and the tests, and run `web/deploy.sh` again.
