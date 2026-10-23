# Installing the Codon Optimiser on your own server

For self-hosting on a Linux server (Ubuntu/Debian assumed). For Streamlit
Community Cloud, see `DEPLOYMENT.md` instead.

If you are Claude Code installing this: work through the steps in order,
run the check at the end of each one, and stop and ask the user before
anything needing sudo, opening firewall ports, or editing an existing nginx
site.

## What is in this folder

```
app.py                  Streamlit front end (entry point)
codon_optimizer.py      the bioinformatics functions
theme.py                palette, CSS, decorative SVG
requirements.txt        streamlit, pandas
.streamlit/config.toml  theme + headless mode - keep it, or the app renders white
README.md, DEPLOYMENT.md, INSTALL.md
```

No database, no secrets, no API keys. The only outbound network call is the
optional UniProt lookup (`https://rest.uniprot.org`); everything else works
offline.

Tested with Python 3.11, streamlit 1.66.1, pandas 2.2.3 and 3.1.0. Any Python >= 3.9
should work.

---

## Step 1 - Unpack and create a virtual environment

```bash
sudo mkdir -p /opt/codon-optimiser
sudo chown "$USER" /opt/codon-optimiser
unzip codon-optimiser.zip -d /opt/codon-optimiser
cd /opt/codon-optimiser

python3 -m venv .venv
.venv/bin/pip install --upgrade pip
.venv/bin/pip install -r requirements.txt
```

Check: `.venv/bin/python -c "import streamlit, pandas; print('ok')"`

## Step 2 - Run it once by hand

```bash
.venv/bin/python -m streamlit run app.py --server.port 8501 --server.address 127.0.0.1
```

Check, from a second shell: `curl -s http://127.0.0.1:8501/_stcore/health`
should print `ok`. Stop it with Ctrl+C.

Use `--server.address 0.0.0.0` instead only if you want to reach it directly
on port 8501 without the reverse proxy in Step 4.

## Step 3 - Keep it running with systemd

Create `/etc/systemd/system/codon-optimiser.service` (replace `www-data` with
the user that should own the process, and give that user read access to
`/opt/codon-optimiser`):

```ini
[Unit]
Description=Codon Optimiser (Streamlit)
After=network.target

[Service]
User=www-data
WorkingDirectory=/opt/codon-optimiser
ExecStart=/opt/codon-optimiser/.venv/bin/python -m streamlit run app.py --server.port 8501 --server.address 127.0.0.1
Restart=on-failure

[Install]
WantedBy=multi-user.target
```

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now codon-optimiser
```

Check: `systemctl status codon-optimiser` shows `active (running)`, and the
curl health check from Step 2 prints `ok`.

## Step 4 - Put nginx in front (optional, recommended)

Streamlit needs **WebSocket** forwarding. Without the `Upgrade`/`Connection`
headers below, the page loads but stays on "Please wait..." forever.

```nginx
server {
    listen 80;
    server_name codon.example.org;   # your domain

    location / {
        proxy_pass http://127.0.0.1:8501;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_read_timeout 86400;
    }
}
```

```bash
sudo nginx -t && sudo systemctl reload nginx
```

For HTTPS, run `sudo certbot --nginx -d codon.example.org` afterwards.

**Serving under a sub-path** (e.g. `example.org/codon/`): add
`--server.baseUrlPath codon` to the `ExecStart` line and use
`location /codon/` with `proxy_pass http://127.0.0.1:8501/codon/;`.

## Updating

Replace the `.py` files (and `requirements.txt` if it changed), then:

```bash
/opt/codon-optimiser/.venv/bin/pip install -r /opt/codon-optimiser/requirements.txt
sudo systemctl restart codon-optimiser
```

## Troubleshooting

| Symptom | Cause |
|---|---|
| Page stuck on "Please wait..." | nginx not forwarding WebSockets (Step 4 headers) |
| App renders plain white | `.streamlit/config.toml` missing - zip tools often skip dot-folders |
| `ModuleNotFoundError: theme` | `theme.py` missing, or service `WorkingDirectory` wrong |
| UniProt lookup says it cannot reach UniProt | server has no outbound HTTPS; the rest of the app still works |
| `journalctl -u codon-optimiser -f` | live logs for anything else |
