# Deploying the Codon Optimiser to Streamlit Community Cloud

## Files this app needs

```
codon-optimiser/
├── app.py                   # Streamlit front end
├── codon_optimizer.py       # the bioinformatics functions
├── theme.py                 # palette, CSS, decorative SVG
├── requirements.txt         # streamlit, pandas
├── .gitignore               # keeps the journal PDFs out of the repo
├── .streamlit/
│   └── config.toml          # base theme - WITHOUT THIS the app deploys white
└── README.md                # optional but good practice
```

That is the whole app. No database, no secrets, no API keys.

> **Warning.** The working folder also contains ~82 MB of copyrighted journal
> PDFs and course documents. They must not go into a public repository.
> The `.gitignore` is deny-by-default (ignore `*`, then re-allow the app
> files) so `git add .` cannot pick them up by accident. Always check
> `git status` before the first commit anyway.

---

## Step 1 - Run it locally first

```bash
python -m streamlit run app.py
```

Use `python -m streamlit`, not plain `streamlit` - on this machine pip
installed the launcher into a folder that is not on PATH.

It opens http://localhost:8501 in your browser. Check:

- both example buttons load a sequence
- the organism dropdown lists all four hosts, and switching one changes the
  GC window, restriction panel and motif count shown in the sidebar
- picking a bacterial host makes the Kozak check read "not applicable"
- "Optimise" produces a sequence and all four tabs render
- pasting nonsense (e.g. `MKXVZ`) shows a red error, not a crash

Stop it with Ctrl+C.

---

## Step 2 - Put it on GitHub

Streamlit Community Cloud deploys **from a GitHub repository**. There is no
file-upload option.

Create the empty repository on github.com first (no README, no .gitignore -
you are pushing those yourself). Then, in PowerShell:

```powershell
cd C:\Users\student\Desktop\claude

git init
git add .
git status          # <- CHECK: only the 6 app files, no PDFs
```

Only once `git status` looks right:

```powershell
git commit -m "Codon optimiser web app"
git branch -M main
git remote add origin https://github.com/<your-username>/codon-optimiser.git
git push -u origin main
```

On the very first commit git may ask who you are:

```powershell
git config --global user.name  "Your Name"
git config --global user.email "you@example.com"
```

The `git push` opens a browser window to sign in to GitHub (Git Credential
Manager ships with Git for Windows). That is expected - no password is typed
into the terminal.

**Public or private?** The free tier allows both, but public is simpler and
has no app limit. Nothing in this code is secret.

---

## Step 3 - Deploy

1. Go to https://share.streamlit.io
2. Sign in **with GitHub** and authorise Streamlit to read your repos
3. Click **New app** -> **Deploy a public app from GitHub**
4. Fill in:
   - Repository: `<your-username>/codon-optimiser`
   - Branch: `main`
   - Main file path: `app.py`
5. Click **Deploy**

First build takes 2-5 minutes while it installs streamlit and pandas.
You get a URL like `https://<your-app-name>.streamlit.app`.

---

## Step 4 - After it is live

**Updating:** push to `main` and the app redeploys automatically. No redeploy
button needed.

**If the build fails:** click "Manage app" (bottom right) for the logs. The
usual causes are a typo in `requirements.txt` or a module that is not listed
there.

**Sleeping:** free apps sleep after about a week with no visitors. The first
visitor after that waits ~30 seconds while it wakes. This is normal and not a
failure - worth knowing before a live demo.

**Resource limits:** 1 GB RAM on the free tier. This app uses a few MB; the
3000 aa cap in `app.py` keeps worst-case runtime around 4 seconds.

---

## Before you share the link

- [ ] Verify the example protein sequences against UniProt (P42212, P01308).
      They were typed from memory and have not been checked against the
      database.
- [ ] Read the LICENSING NOTES at the bottom of `codon_optimizer.py`. The
      Kazusa attribution is already in the app footer; keep it there.
- [ ] Add a one-line disclaimer if anyone outside the course might use it:
      this is a teaching tool, not validated for constructs going to
      synthesis.

---

## Optional: a .gitignore

```
__pycache__/
*.pyc
.venv/
.streamlit/secrets.toml
```
