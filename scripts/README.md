# Scripts

Utility scripts for muset development and maintenance.

## `update_bioconda.sh`

Automates the process of updating the muset package on Bioconda.

### Usage

```bash
# Update to latest git tag
./scripts/update_bioconda.sh

# Or specify a version
./scripts/update_bioconda.sh 0.6.0
```

### Requirements

- Git
- GitHub CLI (`gh`) - optional but recommended for automatic PR creation
- Bioconda-recipes fork at `~/bioconda-recipes`

### First-time setup

1. Fork https://github.com/bioconda/bioconda-recipes
2. Clone your fork: `git clone https://github.com/YOUR_USERNAME/bioconda-recipes.git ~/bioconda-recipes`
3. Add upstream: `cd ~/bioconda-recipes && git remote add upstream https://github.com/bioconda/bioconda-recipes.git`
4. Install GitHub CLI: `conda install gh --channel conda-forge`
5. Authenticate: `gh auth login`

### What it does

1. Syncs your fork with upstream bioconda-recipes
2. Creates a new branch for the update
3. Updates `meta.yaml` with the new version
4. Commits and pushes changes to your fork
5. Creates a pull request to bioconda/bioconda-recipes (if `gh` is installed)

## Automated Updates

The repository also includes a GitHub Actions workflow (`.github/workflows/update-bioconda.yml`) that automatically creates Bioconda update PRs when you publish a GitHub release.

See the main README for more details on automation.
