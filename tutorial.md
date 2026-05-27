# Tutorial: Building a CI/CD Pipeline with GitHub Actions

## What you will build

A GitHub Actions pipeline that runs automatically on every push to `main` — but only when files outside `tests/` are changed. It will cover five stages: linting, automated tests, SBOM and vulnerability scanning, static security analysis, and publishing a binary.

Work through each step, push, and verify the pipeline turns green before moving on.

---

## Setup

GitHub Actions reads workflow files from `.github/workflows/`.  
Create `.github/workflows/pipeline.yml` and start with this skeleton — it already handles the trigger and path filter:

```yaml
name: CI/CD Pipeline

on:
  push:
    branches: [main]
    paths-ignore:
      - 'tests/**'

jobs:
  # your jobs go here
```

`paths-ignore` means the pipeline is skipped entirely when only test files changed.

Push this file. An entry should appear in the **Actions** tab (doing nothing yet — that's fine).

---

## Step 1 — Linting

**What & why:** A linter reads your code without running it and flags style violations and common mistakes. Running it in CI means bad style never silently enters the main branch.

**Tool:** [ruff](https://docs.astral.sh/ruff/) — a modern Python linter and formatter that is significantly faster than flake8 and covers the same rules. Install it with `pip install ruff` and run it with `ruff check .`.

**Your task:** Add a job named `lint` that:
1. checks out the code
2. sets up Python 3.10
3. installs and runs ruff

> The existing codebase uses star imports and long lines. You will need a small `ruff.toml` config file at the project root to keep the linter happy — check the ruff docs for `ignore` and `line-length`.

Push and verify the job is green.

---

## Step 2 — Automated Tests

**What & why:** Running the test suite in CI ensures that every pushed commit is verified against the full test suite, not just the developer's local machine.

**Tool:** [pytest](https://docs.pytest.org/) — already configured in `pytest.ini`.

**Your task:** Add a job named `test` that installs dependencies from `requirements.txt` and runs `pytest tests/ -v`. It should only start after `lint` has passed — use the `needs:` key.

> One CLI command in this project depends on **Graphviz**, a system-level binary not installable via pip. Without it, three tests are silently skipped. To run those too, install Graphviz with `apt-get` before the pip step.

Push and verify both `lint` and `test` are green.

---

## Step 3 — SBOM & Vulnerability Scan

**What & why:** An SBOM (Software Bill of Materials) is a machine-readable inventory of every library your software depends on. Scanning it against a CVE database lets you catch known vulnerabilities in your dependencies before they reach production.

**Tools:**
- [anchore/sbom-action](https://github.com/anchore/sbom-action) — generates the SBOM (no extra tooling needed, just a `uses:` line)
- [anchore/scan-action](https://github.com/anchore/scan-action) — scans the SBOM for CVEs using [Grype](https://github.com/anchore/grype)

**Your task:** Add a job named `sbom` that generates a CycloneDX JSON SBOM, scans it, and uploads it as a pipeline artifact so it can be downloaded from the Actions UI. It should only start after `test` has passed.

> The dependencies in `requirements.txt` are pinned to 2022 versions and likely have known CVEs. Neither the SBOM scan nor the pip-audit scan should block the pipeline — use `fail-build: false` and `continue-on-error: true` so findings are visible in the log without stopping the build. This makes for a great live discussion about dependency hygiene.

Push, open **Actions → your run → Artifacts**, and download the generated SBOM to see what it contains.

---

## Step 4 — Static Security Scan (SAST)

**What & why:** Unlike a vulnerability scanner (which checks your *dependencies*), a SAST tool reads *your own source code* and flags dangerous patterns — hard-coded secrets, use of unsafe functions, injection risks, and more.

**Tool:** [bandit](https://bandit.readthedocs.io/) — the de-facto Python SAST tool. Install it with `pip install bandit` and run it with `bandit -r . --exclude .venv,tests`.

**Your task:** Add a job named `security` that installs and runs bandit. Use a flag to set a minimum severity threshold so the job fails only on findings above a certain level. Like `sbom`, it should only start after `test` has passed.

Push and read the bandit output — it prints a readable summary of any findings.

> **Alternative:** [Semgrep](https://semgrep.dev/docs/semgrep-ci/sample-ci-configs/#github-actions) has a ready-made GitHub Action and supports many languages beyond Python.

---

## Step 5 — Build & Publish a Binary

**What & why:** The final stage packages the tool into a standalone executable that anyone can run without installing Python. This is the "ship it" step — it only runs once everything else has passed.

**Tool:** [PyInstaller](https://pyinstaller.org/en/stable/) — bundles a Python script and all its dependencies into a single binary. Install it with `pip install pyinstaller` and run it with `pyinstaller --onefile main.py --name optimizhelper`.

**Your task:** Add a job named `publish` that:
1. builds the binary
2. uploads it as a pipeline artifact named `optimizhelper-linux`
3. only starts **after** `lint`, `test`, `sbom`, and `security` have all passed — use the [`needs:`](https://docs.github.com/en/actions/writing-workflows/workflow-syntax-for-github-actions#jobsjob_idneeds) key for this

Push, then go to **Actions → your run → Artifacts** and download the binary. Try running it.

---

---

## Possible Solutions

Stuck? Expand the section for the step you are on.

<details>
<summary><strong>Step 1 — Linting</strong></summary>

Add `ruff.toml` at the project root:

```toml
line-length = 120
exclude = [".venv", "tests"]

[lint]
ignore = ["E402", "F401", "F403", "F405"]
```

Job:

```yaml
  lint:
    name: Lint
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - run: pip install ruff
      - run: ruff check .
```

</details>

<details>
<summary><strong>Step 2 — Tests</strong></summary>

```yaml
  test:
    name: Test
    runs-on: ubuntu-latest
    needs: [lint]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - name: Install system dependencies
        run: sudo apt-get install -y graphviz
      - name: Install Python dependencies
        run: pip install -r requirements.txt pytest
      - name: Run test suite
        run: pytest tests/ -v
```

</details>

<details>
<summary><strong>Step 3 — SBOM & Vulnerability Scan</strong></summary>

```yaml
  sbom:
    name: SBOM & Vulnerability Scan
    runs-on: ubuntu-latest
    needs: [test]
    steps:
      - uses: actions/checkout@v4

      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'

      - name: Install Syft
        run: |
          curl -sSfL https://raw.githubusercontent.com/anchore/syft/main/install.sh | sh -s -- -b /usr/local/bin

      - name: Generate SBOM
        run: syft . -o cyclonedx-json=sbom.json

      - name: Upload SBOM as artifact
        uses: actions/upload-artifact@v4
        with:
          name: sbom
          path: sbom.json

      - name: Scan SBOM for vulnerabilities
        uses: anchore/scan-action@v3
        with:
          sbom: sbom.json
          fail-build: false
        continue-on-error: true

      - name: Scan dependencies for vulnerabilities
        run: |
          pip install pip-audit
          pip-audit -r requirements.txt
        continue-on-error: true
```

`anchore/scan-action` attempts a CVE scan of the SBOM but won't fail the build — it can be brittle due to database caching issues, so `continue-on-error: true` is a safety net. `pip-audit` then does a reliable scan directly against `requirements.txt`. Both steps report findings in the log without blocking the pipeline — the old pinned dependencies have known CVEs, which makes for a great live discussion about dependency hygiene.

</details>

<details>
<summary><strong>Step 4 — Static Security Scan</strong></summary>

```yaml
  security:
    name: Static Security Scan
    runs-on: ubuntu-latest
    needs: [test]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - run: pip install bandit
      - name: Run bandit
        run: bandit -r . --exclude .venv,tests -ll
```

`-ll` reports medium severity and above. Use `-l` to also include low-severity findings.

</details>

<details>
<summary><strong>Step 5 — Build & Publish</strong></summary>

```yaml
  publish:
    name: Build & Publish Binary
    runs-on: ubuntu-latest
    needs: [lint, test, sbom, security]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - name: Install dependencies
        run: pip install -r requirements.txt pyinstaller
      - name: Build standalone binary
        run: pyinstaller --onefile main.py --name optimizhelper
      - name: Upload binary as artifact
        uses: actions/upload-artifact@v4
        with:
          name: optimizhelper-linux
          path: dist/optimizhelper
```

</details>

<details>
<summary><strong>Complete <code>pipeline.yml</code></strong></summary>

```yaml
name: CI/CD Pipeline

on:
  push:
    branches: [main]
    paths-ignore:
      - 'tests/**'

jobs:

  lint:
    name: Lint
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - run: pip install ruff
      - run: ruff check .

  test:
    name: Test
    runs-on: ubuntu-latest
    needs: [lint]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - name: Install system dependencies
        run: sudo apt-get install -y graphviz
      - name: Install Python dependencies
        run: pip install -r requirements.txt pytest
      - name: Run test suite
        run: pytest tests/ -v

  sbom:
    name: SBOM & Vulnerability Scan
    runs-on: ubuntu-latest
    needs: [test]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - name: Generate SBOM
        uses: anchore/sbom-action@v0
        with:
          format: cyclonedx-json
      - name: Upload SBOM as artifact
        uses: actions/upload-artifact@v4
        with:
          name: sbom
          path: OptimizHelper-sbom.cyclonedx.json
      - name: Scan SBOM for vulnerabilities
        uses: anchore/scan-action@v3
        with:
          sbom: OptimizHelper-sbom.cyclonedx.json
          fail-build: false
        continue-on-error: true
      - name: Scan dependencies for vulnerabilities
        run: |
          pip install pip-audit
          pip-audit -r requirements.txt
        continue-on-error: true

  security:
    name: Static Security Scan
    runs-on: ubuntu-latest
    needs: [test]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - run: pip install bandit
      - name: Run bandit
        run: bandit -r . --exclude .venv,tests -ll

  publish:
    name: Build & Publish Binary
    runs-on: ubuntu-latest
    needs: [lint, test, sbom, security]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - name: Install dependencies
        run: pip install -r requirements.txt pyinstaller
      - name: Build standalone binary
        run: pyinstaller --onefile main.py --name optimizhelper
      - name: Upload binary as artifact
        uses: actions/upload-artifact@v4
        with:
          name: optimizhelper-linux
          path: dist/optimizhelper
```

</details>
