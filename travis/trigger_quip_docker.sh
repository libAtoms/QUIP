#!/bin/bash

curl -s --output /dev/null -X POST \
    -H "Content-Type: application/json" \
    -H "Accept: application/json" \
    -H "Travis-API-Version: 3" \
    -H "Authorization: token ${TRAVIS_TOKEN}" \
    -d "{\"request\": {\"branch\": \"master\"}}" \
    https://api.travis-ci.org/repo/libatoms%2Fquip-docker/requests