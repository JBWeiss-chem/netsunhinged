# NetsUnhinged - production image (any VPS, university VM or container platform)
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PORT=8765

# unprivileged user
RUN useradd --create-home --uid 1000 app

# install the package itself (front end and examples are package data), then drop the sources
COPY pyproject.toml README.md LICENSE /tmp/src/
COPY netsunhinged /tmp/src/netsunhinged
RUN pip install --no-cache-dir "/tmp/src[server]" && rm -rf /tmp/src

USER app
WORKDIR /home/app
EXPOSE 8765

# exactly ONE worker process: sessions and analysis jobs are kept in memory;
# threads serve concurrent requests, analyses run in the app's own worker pool.
# No access log (it would record visitors' IP addresses); errors still go to stderr.
CMD gunicorn --workers 1 --threads 8 --timeout 120 --bind 0.0.0.0:${PORT} \
    netsunhinged.server:app
