#!/bin/bash

echo $github_deploy_key | base64 -d >~/.ssh/github_deploy_key
chmod 0600 ~/.ssh/github_deploy_key
eval $(ssh-agent -s)
ssh-add ~/.ssh/github_deploy_key
git clone --depth 10 git@github.com:libAtoms/GAP.git src/GAP
