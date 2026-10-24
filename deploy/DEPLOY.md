# Deploying NetsUnhinged

The app is one Python process (gunicorn, **one worker**, several threads). Sessions and analysis jobs
are kept in memory, so never run more than one worker or replica.

## Before you publish

* `papers/` (copyrighted articles) and `test_cifs/` (unpublished structures) are excluded by
  `.dockerignore`; check this if you deploy in another way.
* The example CIFs from the CSD carry the CCDC data-access notice. Check the CCDC terms before
  redistributing them publicly, or replace them with links. The IZA zeolite files are freely available.
* Depending on where you operate the site, the law may require a legal notice or a privacy notice
  (e.g. an *Impressum* and *Datenschutzerklärung* in Germany). Check the rules that apply to you.
* What the app does with personal data (for a privacy notice):
  * Uploaded structures are kept in memory only, are deleted after 1 h of inactivity (earlier when
    `NU_MAX_SESSIONS` is reached), and are never written to disk or logged.
  * The app writes no access log: the Docker image runs gunicorn without `--access-logfile`. Internal
    errors are logged with a reference number, the request path and a traceback, without IP addresses or
    file contents.
  * For rate limiting, client IP addresses are held in memory for at most 60 s.
  * No cookies, no browser storage, no third-party requests: three.js is served from
    `netsunhinged/static/vendor/three/` (version 0.160.0, MIT licence), not from a CDN.
  * The hosting provider or reverse proxy may keep its own connection logs (Caddy only if you enable
    its `log` directive). Mention the hoster in the privacy notice.

## Limits (environment variables)

| Variable | Default | Meaning |
|---|---|---|
| `NU_MAX_UPLOAD_MB` | 5 | maximum request size |
| `NU_MAX_ATOMS` | 6000 | maximum atoms per unit cell after symmetry expansion |
| `NU_MAX_DOF` | 6000 | larger models skip the infinitesimal analysis |
| `NU_WORKERS` | 1 | analyses running in parallel |
| `NU_MAX_QUEUE` | 4 | analyses waiting; further requests get "server busy" |
| `NU_JOB_TIMEOUT` | 600 | seconds per analysis |
| `NU_HEAVY_CONCURRENCY` | 2 | parallel load/model/export requests |
| `NU_TTL` | 3600 | idle seconds before a session is deleted |
| `NU_MAX_SESSIONS` | 200 | stored structures / jobs (oldest dropped first) |
| `NU_RATE_PER_MIN` | 60 | heavy requests (load, model, analyse, export) per client and minute |
| `NU_RATE_FRAMES_PER_MIN` | 600 | deformed-frame requests (regime slider) per client and minute |
| `NU_TRUST_PROXY` | 0 | set to 1 behind a reverse proxy, so rate limits use `X-Forwarded-For` |
| `NU_FRAME_ANCESTORS` | `'self'` | who may embed the page in a frame (CSP `frame-ancestors`) |

## Your own server (VPS or university VM)

Use the published image (built by GitHub Actions for every push to `main` and every release)…

```bash
docker run -d --name netsunhinged --restart unless-stopped \
  -p 127.0.0.1:8765:8765 -e NU_TRUST_PROXY=1 \
  ghcr.io/jbweiss-chem/netsunhinged:latest
```

…or build it yourself from a checkout: `docker build -t netsunhinged .`, then the same command with
`netsunhinged` as the image name.

The container listens only on localhost. Put Caddy in front: it gets and renews the HTTPS
certificate automatically. `deploy/Caddyfile`:

```
netsunhinged.example.org {
    encode gzip
    reverse_proxy 127.0.0.1:8765
}
```

Firewall: open only ports 22, 80 and 443. Keep the OS and the image up to date (rebuild regularly).

Without Docker: `pip install "netsunhinged[server]"` (or `pip install ".[server]"` in a checkout), then run
`gunicorn --workers 1 --threads 8 --bind 127.0.0.1:8765 netsunhinged.server:app` as an unprivileged
user (e.g. via systemd), behind the same Caddy configuration.
