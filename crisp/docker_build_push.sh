#!/usr/bin/env bash
if docker build -t ah-crisp-test .  ; then
  docker tag ah-crisp-test gcr.io/fdl-us-astronaut-health/ah-crisp-test
  docker push gcr.io/fdl-us-astronaut-health/ah-crisp-test

#  docker tag ah-causal-ensemble-docker registry.gitlab.com/frontierdevelopmentlab/astronaut-health/ah-causal-ensemble
#  docker push registry.gitlab.com/frontierdevelopmentlab/astronaut-health/ah-causal-ensemble
fi
