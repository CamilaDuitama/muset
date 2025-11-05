#!/bin/bash
# Script to automate Bioconda recipe updates for muset
# Usage: ./update_bioconda.sh [version]

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get version from argument or latest git tag
if [ -z "$1" ]; then
    VERSION=$(git describe --tags --abbrev=0 | sed 's/^v//')
    echo -e "${YELLOW}No version specified, using latest tag: $VERSION${NC}"
else
    VERSION=$1
    echo -e "${GREEN}Using specified version: $VERSION${NC}"
fi

# Configuration
GITHUB_USER="CamilaDuitama"
BIOCONDA_FORK="$HOME/bioconda-recipes"
BRANCH_NAME="muset-${VERSION}"

echo -e "\n${GREEN}=== Bioconda Recipe Update Script ===${NC}\n"

# Step 1: Check if bioconda-recipes fork exists
if [ ! -d "$BIOCONDA_FORK" ]; then
    echo -e "${YELLOW}Bioconda recipes not found. Cloning your fork...${NC}"
    cd $(dirname "$BIOCONDA_FORK")
    git clone "https://github.com/${GITHUB_USER}/bioconda-recipes.git"
    cd "$BIOCONDA_FORK"
    git remote add upstream https://github.com/bioconda/bioconda-recipes.git
else
    cd "$BIOCONDA_FORK"
    echo -e "${GREEN}✓ Found bioconda-recipes at $BIOCONDA_FORK${NC}"
fi

# Step 2: Update fork
echo -e "\n${YELLOW}Updating fork with upstream changes...${NC}"
git fetch upstream
git checkout master
git merge upstream/master --ff-only || {
    echo -e "${RED}Failed to update master. Please resolve conflicts manually.${NC}"
    exit 1
}
git push origin master

# Step 3: Create new branch
echo -e "\n${YELLOW}Creating branch: $BRANCH_NAME${NC}"
git checkout -b "$BRANCH_NAME" 2>/dev/null || git checkout "$BRANCH_NAME"

# Step 4: Update meta.yaml
echo -e "\n${YELLOW}Updating recipes/muset/meta.yaml...${NC}"
RECIPE_FILE="recipes/muset/meta.yaml"

if [ ! -f "$RECIPE_FILE" ]; then
    echo -e "${RED}Error: $RECIPE_FILE not found!${NC}"
    exit 1
fi

# Backup original file
cp "$RECIPE_FILE" "${RECIPE_FILE}.backup"

# Update version
sed -i "s/{% set version = \".*\" %}/{% set version = \"${VERSION}\" %}/" "$RECIPE_FILE"

# Reset build number to 0 for new version
sed -i "s/  number: .*/  number: 0/" "$RECIPE_FILE"

# Show changes
echo -e "\n${GREEN}Changes made:${NC}"
diff "${RECIPE_FILE}.backup" "$RECIPE_FILE" || true
rm "${RECIPE_FILE}.backup"

# Step 5: Commit changes
echo -e "\n${YELLOW}Committing changes...${NC}"
git add "$RECIPE_FILE"
git commit -m "Update muset to ${VERSION}" || {
    echo -e "${YELLOW}No changes to commit (recipe already up to date?)${NC}"
    exit 0
}

# Step 6: Push to fork
echo -e "\n${YELLOW}Pushing to your fork...${NC}"
git push origin "$BRANCH_NAME" || {
    echo -e "${RED}Failed to push. Do you have write access to your fork?${NC}"
    exit 1
}

# Step 7: Create PR (requires GitHub CLI)
echo -e "\n${YELLOW}Creating Pull Request...${NC}"
if command -v gh &> /dev/null; then
    gh pr create \
        --repo bioconda/bioconda-recipes \
        --head "${GITHUB_USER}:${BRANCH_NAME}" \
        --title "Update muset to ${VERSION}" \
        --body "This PR updates \`muset\` to version ${VERSION}.

## Changes in v${VERSION}
See release notes: https://github.com/${GITHUB_USER}/muset/releases/tag/v${VERSION}

## Checklist
- [x] Version updated in meta.yaml
- [x] Build number reset to 0
- [ ] CI tests passed (automated)

---
*This PR was created using the update_bioconda.sh script.*"
    
    echo -e "\n${GREEN}✓ Pull Request created successfully!${NC}"
    echo -e "${GREEN}View it at: https://github.com/bioconda/bioconda-recipes/pulls${NC}"
else
    echo -e "${YELLOW}GitHub CLI (gh) not found. Please create PR manually:${NC}"
    echo -e "${GREEN}1. Go to: https://github.com/${GITHUB_USER}/bioconda-recipes${NC}"
    echo -e "${GREEN}2. Click 'Compare & pull request' for branch: ${BRANCH_NAME}${NC}"
    echo -e "${GREEN}3. Create PR to: bioconda/bioconda-recipes:master${NC}"
fi

echo -e "\n${GREEN}=== Done! ===${NC}\n"
