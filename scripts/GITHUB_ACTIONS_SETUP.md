# Setting Up Bioconda Auto-Update

## Steps to Enable GitHub Actions Automation

### 1. Create a Personal Access Token (PAT)

1. Go to: https://github.com/settings/tokens/new
2. Settings:
   - **Note:** `BIOCONDA_PAT`
   - **Expiration:** No expiration (or 1 year)
   - **Scopes:** Check these boxes:
     - ✅ `repo` (Full control of private repositories)
     - ✅ `workflow` (Update GitHub Action workflows)
3. Click **"Generate token"**
4. **Copy the token** (you won't see it again!)

### 2. Add Token to Repository Secrets

1. Go to: https://github.com/CamilaDuitama/muset/settings/secrets/actions
2. Click **"New repository secret"**
3. Settings:
   - **Name:** `BIOCONDA_PAT`
   - **Secret:** (paste your token from step 1)
4. Click **"Add secret"**

### 3. Fork bioconda-recipes (if not done)

1. Go to: https://github.com/bioconda/bioconda-recipes
2. Click **"Fork"** (top right)
3. This creates: https://github.com/CamilaDuitama/bioconda-recipes

### 4. Test the Workflow

After creating a release, check:
- https://github.com/CamilaDuitama/muset/actions
- Look for "Update Bioconda Recipe" workflow
- If it fails, check the logs for errors

## Troubleshooting

### Workflow doesn't run
- Check that `BIOCONDA_PAT` secret exists in repo settings
- Verify the workflow file is in `.github/workflows/update-bioconda.yml`
- Check that you published a release (not just pushed a tag)

### PR creation fails
- Ensure you forked bioconda-recipes
- Verify PAT has `repo` and `workflow` scopes
- Check PAT hasn't expired

### Can't find the PR
- Search: https://github.com/bioconda/bioconda-recipes/pulls?q=muset
- Check your fork: https://github.com/CamilaDuitama/bioconda-recipes/branches

## Manual Alternative

If automation fails, use the manual script:
```bash
./scripts/update_bioconda.sh 0.6.0
```
